function OUTPUT=EIOT_exclusive_signature_PB_star1(S,S_I,SPEC_test,SPEC_target)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code will search the number of best residual spectra to reach less 
% than 1% change of sse for both EIOT and PACLS.
%
% SPEC_CAL, S SPEC_sim and SPEC_test need to be preprocessed so that we can
% switch back and forth between snv and derivative
%
% SPEC_test has to be one spectrum a time
%
% Given project-dependent S_I is better than project-independent ones, let
% us always use project dependent S_I here.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EIOT
disp('EIOT section begins...')
[c_E_hat,Em,sse_old,lamda] = EIOT_PRED([S'],0,snv(SPEC_test)',[-inf inf],1);    %1st

old=S;
new=S_I;

for j=1:size(S_I,1)
    [c_E_hat,Em,sse_new(j),lamda] = EIOT_PRED([old' new(j,:)'],size([old;new(j,:)],1)-size(S,1),(SPEC_test)',[-inf inf],1);    %2nd
end
[C,I]=sort(sse_new,'ascend');


vec=new(I(1),:);
old=[old;new(I(1),:)];
new=new(I(2:end),:);

while (min(sse_old)-min(sse_new))/min(sse_old)*100>1
    sse_old=sse_new;
    sse_new=[];
    for j=1:size(new,1)
        [c_E_hat,Em,sse_new(j),lamda] = EIOT_PRED([old' new(j,:)'],size([old;new(j,:)],1)-size(S,1),(SPEC_test)',[-inf inf],1);    %3rd
    end

    [C,I]=sort(sse_new,'ascend');

    % sse_new
    if isempty(I)
        break
    else
        old=[old;new(I(1),:)];
        vec(size(vec,1)+1,:)=new(I(1),:);
        new=new(I(2:end),:);
    end
  
end
% (min(sse_old)-min(sse_new))/min(sse_old)*100


[c_E_hat,Em,sse_new,lamda] = EIOT_PRED([S' vec(1:end-1,:)'],size(vec,1)-1,(SPEC_test)',[-inf inf],1);

OUTPUT.r_hat=c_E_hat';
OUTPUT.EIOT_res=Em';
OUTPUT.lamda=lamda;
sse_new
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PACLS
disp('PACLS section begins...')
Y_PRED_PACLS=(SPEC_test)*pinv([S]);   %1st
RES=Y_PRED_PACLS*[S]-SPEC_test;
sse_old=sum(RES.^2);

old=S;
new=S_I;

for j=1:size(S_I,1)
    Y_PRED_PACLS=(SPEC_test)*pinv([old;new(j,:)]);    %2nd
    RES=Y_PRED_PACLS*[old;new(j,:)]-SPEC_test;
    sse_new(j)=sum(RES.^2);
end
[C,I]=sort(sse_new,'ascend');

old=[old;new(I(1),:)];
vec_pacls=new(I(1),:);
new=new(I(2:end),:);

while (min(sse_old)-min(sse_new))/min(sse_old)*100>1
    sse_old=sse_new;
    sse_new=[];
    for j=1:size(new,1)
        Y_PRED_PACLS=(SPEC_test)*pinv([old;new(j,:)]);   %3rd
        RES=Y_PRED_PACLS*[old;new(j,:)]-SPEC_test;
        sse_new(j)=sum(RES.^2);
    end
    [C,I]=sort(sse_new,'ascend');
    % sse_new
    if isempty(I)
        break
    else
        old=[old;new(I(1),:)];
        vec_pacls(size(vec_pacls,1)+1,:)=new(I(1),:);
        new=new(I(2:end),:);
    end

end
% (min(sse_old)-min(sse_new))/min(sse_old)*100


Y_PRED_PACLS=(SPEC_test)*pinv([S;vec_pacls(1:end-1,:)]);
RES=Y_PRED_PACLS*[S;vec_pacls(1:end-1,:)]-SPEC_test;

OUTPUT.PACLS_pred=Y_PRED_PACLS;
OUTPUT.PACLS_res=RES;
sum(RES.^2)

OUTPUT_ideal = NAS_FILTER_1([],(S(2:end,:)),(S(1,:)),[],size(S,1)-1);

SSI=sum(SPEC_test.^2);
SSR_EIOT=sum(OUTPUT.EIOT_res.^2);
SSR_PACLS=sum(OUTPUT.PACLS_res.^2);
OUTPUT_mix = NAS_FILTER_1([],(S(2:end,:)),SPEC_test,[],size(S,1)-1);
temp=(OUTPUT_mix.X_PREP{1,end}'*OUTPUT_ideal.X_PREP{1,end})/(norm(OUTPUT_mix.X_PREP{1,end}')*norm(OUTPUT_ideal.X_PREP{1,end}'));
OUTPUT.NAS_T=acosd(temp);

OUTPUT_target = NAS_FILTER_1([],(S(2:end,:)),SPEC_target,[],size(S,1)-1);

temp1=abs((norm(OUTPUT_mix.X_PREP{1,end}')./norm(OUTPUT_ideal.X_PREP{1,end}'))-(norm(OUTPUT_target.X_PREP{1,end}')./norm(OUTPUT_ideal.X_PREP{1,end}')));
OUTPUT.SND=temp1;

OUTPUT.EIOT_ssr=SSR_EIOT;
OUTPUT.PACLS_ssr=SSR_PACLS;
OUTPUT.EIOT_ssr_relative=SSR_EIOT./SSI*100;
OUTPUT.PACLS_ssr_relative=SSR_PACLS./SSI*100;
OUTPUT.vec_EIOT=vec;
OUTPUT.vec_PACLS=vec_pacls;

