function OUTPUT=EIOT_exclusive_signature_PB_star3(S,S_I,SPEC_test,SPEC_sim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code will the entirety of independent residual spectrum for both 
% EIOT and PACLS.
%
% S, SPEC_sim and SPEC_test need to be preprocessed so that we can
% switch back and forth between snv and derivative
%
% SPEC_test does not have to be one spectrum a time
%
% Given project-dependent S_I is better than project-independent ones, let
% us always use project dependent S_I here.
%
for i=1:size(SPEC_test,1)

    [c_E_hat,Em,sse,lamda] = EIOT_PRED([S' S_I'],6,SPEC_test(i,:)',[-inf inf],1); %that number needs to be changed per project
    
    OUTPUT.EIOT_pred(i,:)=c_E_hat';
    OUTPUT.EIOT_res(i,:)=Em;
    OUTPUT.EIOT_ssr(i)=sse;
    OUTPUT.lamda_eqlin=lamda.eqlin;
    OUTPUT.lamda_lower(i,:)=lamda.lower(1:size(S,1))';
    OUTPUT.lamda_upper(i,:)=lamda.upper(1:size(S,1))';
end

Y_PRED_PACLS=(SPEC_test)*pinv([S;S_I]);
RES=Y_PRED_PACLS*[S;S_I]-SPEC_test;


OUTPUT.PACLS_pred=Y_PRED_PACLS;
OUTPUT.PACLS_res=RES;

OUTPUT_ideal = NAS_FILTER_1([],(S(2:end,:)),(SPEC_sim),[],size(S,1)-1);

for i=1:size(SPEC_test,1)
    SSI(i)=sum(SPEC_test(i,:).^2);
    SSR_PACLS(i)=sum(OUTPUT.PACLS_res(i,:).^2);
    OUTPUT_mix = NAS_FILTER_1([],(S(2:end,:)),SPEC_test(i,:),[],size(S,1)-1);
    temp=(OUTPUT_mix.X_PREP{1,end}'*OUTPUT_ideal.X_PREP{1,end})/(norm(OUTPUT_mix.X_PREP{1,end}')*norm(OUTPUT_ideal.X_PREP{1,end}'));
    OUTPUT.NAS_T(i)=acosd(temp);
end

OUTPUT.PACLS_ssr=sse;
OUTPUT.EIOT_ssr_relative=OUTPUT.EIOT_ssr./SSI*100;
OUTPUT.PACLS_ssr_relative=SSR_PACLS./SSI*100;


