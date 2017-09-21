%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = modelCorrections(model)
% Corrects various issues in yeast 7.
%
% Benjamín J. Sánchez. Last edited: 2015-12-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = modelCorrections(model)

%Correct glucan coefficients in biomass reaction:
model.S(strcmp(model.mets,'s_0002'),strcmp(model.rxns,'r_4041')) = 0;
model.S(strcmp(model.mets,'s_0001'),strcmp(model.rxns,'r_4041')) = -0.8506;
model.S(strcmp(model.mets,'s_0004'),strcmp(model.rxns,'r_4041')) = -0.2842;

%Correctly represent proton balance inside cell:
model.lb(strcmp(model.rxns,'r_1824')) = 0;  %Block free H+ export
model.ub(strcmp(model.rxns,'r_1250')) = 0;  %Block free putrescine export
model.ub(strcmp(model.rxns,'r_1259')) = 0;  %Block free spermidine export

%CHANGES IN OX.PHO.:
%COMPLEX III: H+ pumping according to Förster 2003 (75% eff in P/O ratio):
model.S(strcmp(model.mets,'s_0799'),strcmp(model.rxns,'r_0439')) = -1.5;
model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_0439')) = +3;
%COMPLEX IV: H+ pumping according to Förster 2003 (75% eff in P/O ratio):
model.S(strcmp(model.mets,'s_0799'),strcmp(model.rxns,'r_0438')) = -6;
model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_0438')) = +3;
%COMPLEX IV: Normalize rxn by the number of ferrocytochromes c:
rxn_pos            = strcmp(model.rxns,'r_0438');
ferro_S            = abs(model.S(strcmp(model.mets,'s_0710'),rxn_pos));
model.S(:,rxn_pos) = model.S(:,rxn_pos)./ferro_S;
%COMPLEX V: For 1 ATP 3 H+ are needed, not 4:
model.S(strcmp(model.mets,'s_0799'),strcmp(model.rxns,'r_0226')) = +2;
model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_0226')) = -3;

%Correct rev vector: true if LB < 0 & UB > 0:
model.rev = boolean((model.lb < 0).*(model.ub > 0));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%