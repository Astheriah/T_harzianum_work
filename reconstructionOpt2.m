function modelo_final=reconstructionOpt2(model_template,modelo_id,blastStructureCode,maxE,minLen,minSim)

%
model_template.id=model_template.description;


models= {model_template};
getModelFor=modelo_id;
preferredOrder={};
strictness=1;
onlyGenesInModels=false();
mapNewGenesToOld=true();

draftModel=getModelFromHomology(models,blastStructureCode,getModelFor,preferredOrder,strictness,onlyGenesInModels,maxE,minLen,minSim,mapNewGenesToOld);
model_draft1=draftModel;
modelo_templado=model_template;

 for i=1:length(model_draft1.rxns) %<- cambiar grRules del draft al modelo que uses como templado
     p=strmatch(model_draft1.rxns(i,1),modelo_templado.rxns);
     modelo_templado.grRules(p,1)=model_draft1.grRules(i,1);
 
 end  
modelo_templado =rmfield(modelo_templado,'rules');
modelo_templado = buildRxnGeneMat(modelo_templado);
modelo_templado = generateRules(modelo_templado);
modelo_templado = removeUnusedGenes(modelo_templado);

%modelo_templado.genes=unique(vertcat(modelo_templado.genes,draftModel.genes));

modelo_optimiza=modelo_templado;
clear Index;
clear Biomassrxn;
ATPindex=strmatch('ATPM', modelo_optimiza.rxns);
if isempty(ATPindex)==1
    rxnName={'ATPM','ATP maintenance requirement'};
            metaboliteList={'atp_c','h2o_c', 'adp_c', 'h_c', 'pi_c'}; %Make sure that the naming is consisten for these metabolites
            stoichCoeffList=[-1,-1,1,1,1];
            revFlag=false;
            lowerBound=0;
            upperBound=1000;
            objCoeff=0;
            subSystem=' ';
            grRule=' ';
            geneNameList='';
            systNameList='';
            checkDuplicate=true;
            confScores=4;
            rxnECNumbers=' ';
            rxnKEGGID=' ';
            rxnReferences=' ';
            rxnNotes=' ';
            modelo_optimiza=addReaction2(modelo_optimiza,rxnName,metaboliteList,stoichCoeffList,revFlag,lowerBound,upperBound,objCoeff,subSystem,grRule,geneNameList,systNameList,checkDuplicate,confScores);
end    
ATPreaction=modelo_optimiza.rxns(ATPindex);
Biomassindex=strfind(modelo_optimiza.rxns,'BIOMASS');
cont3=1;
for i=1:length(Biomassindex)
    if isempty(Biomassindex{i})==0
        Biomassrxn(cont3,1)=modelo_optimiza.rxns(i,1);
        cont3=cont3+1;
    end    

end 
ExCh=findExcRxns(modelo_optimiza);
a=modelo_optimiza.rxns(ExCh);
%match=strfind(modelo_optimiza.subSystems, 'Transport');
cont2=1;
cont=1;
clear Transport_index;
clear transport
% for i=1:length(match)
%     if isempty(match{i})==0
%         Transport_index(cont2,1)=i;
%         transport(cont2,1)=modelo_optimiza.rxns(i,1);
%         transport(cont2,2)=modelo_optimiza.subSystems(i,1);
%         cont2=cont2+1;
% 
%     end    
% 
% end
%TrCh=transport(:,1);
Index=unique(vertcat(Biomassrxn,ATPreaction,model_draft1.rxns,a));
P_Index=setdiff(modelo_templado.rxns,Index);
numel(P_Index)
P_Index=unique(P_Index);
numel(P_Index)
modelo=modelo_optimiza;
ATPindex=strmatch('ATPM', modelo);
if isempty(ATPindex)==1
        rxnName={'ATPM','ATP maintenance requirement'};
        metaboliteList={'atp_c','h2o_c', 'adp_c', 'h_c', 'pi_c'}; %Make sure that the naming is consisten for these metabolites
        stoichCoeffList=[-1,-1,1,1,1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem=' ';
        grRule=' ';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores=4;
        rxnECNumbers=' ';
        rxnKEGGID=' ';
        rxnReferences=' ';
        rxnNotes=' ';
        modelo=addReaction2(modelo,rxnName,metaboliteList,stoichCoeffList,revFlag,lowerBound,upperBound,objCoeff,subSystem,grRule,geneNameList,systNameList,checkDuplicate,confScores);
end  
    
f=strmatch('NTP1',P_Index, 'exact');
if isempty(f)==0
    NTP1(1,1)=P_Index(f,1);
    P_Index=setdiff(P_Index,NTP1);
end    
modelo_remove=modelo;
modeloRep=modelo;
FBAsolution=optimizeCbModel(modeloRep)
FBAsolution.f

    
for i=1:numel(P_Index)   
    
    modelo2=modeloRep;
    %FBAsolution=optimizeCbModel(modeloRep);  
    modelo_remove=removeRxns(modelo2,P_Index(i));

    modelo_remove.csense={};

    for j=1:length(modelo_remove.mets)
    modelo_remove.csense{j}='E';
    end
    modelo_remove.csense=char(modelo_remove.csense);

    if length(modelo_remove.S(:,1))<length(modelo_remove.S(1,:))
       FBAsolution=optimizeCbModel(modelo_remove,'max','one');
       ATPmodelo=modelo_remove;
       ATPmodelo=changeObjective(ATPmodelo,'ATPM');
       Exch3=findExcRxns(ATPmodelo);
       ATPmodelo.lb(Exch3)=0;
       ATPmodelo.ub(Exch3)=1000;

           FBAtrATP=optimizeCbModel(ATPmodelo, 'max', 'one');
%            FBAsolution.f
%            FBAtrATP.f   
           if FBAsolution.f>0.01 && (FBAtrATP.f==0 || isnan(FBAtrATP.f))
               modelo2=modelo_remove;
               modeloRep=modelo_remove;
           else
                modelo2=modeloRep;
               %modeloRep=modeloRep;
           end 
   else
       break
   end
end
%Prueba de calidad de ATP
modelFixed = modelo2;
% cont=1;
% for i=1:length(modelFixed.mets)
%     checkMet=findRxnsFromMets(modelFixed,modelFixed.mets(i));
%     if isempty(checkMet)
%         listDelete{cont,1}=modelFixed.mets{i,1};
%         cont=cont+1;
% 
%     end
% end
%modelFixed = removeMetabolites(modelFixed, listDelete);
modelFixed =rmfield(modelFixed,'rules');
modelFixed = buildRxnGeneMat(modelFixed);
modelFixed = generateRules(modelFixed);
modelFixed = removeUnusedGenes(modelFixed);
modelo_final=modelFixed;
modeloATP=modelo_final;
% modelo_final=buildRxnGeneMat(modelo_final);
Exch2=findExcRxns(modeloATP);
modeloATP=changeObjective(modeloATP,'ATPM');
modeloATP.lb(Exch2)=0;
modeloATP.ub(Exch2)=1000;
ATPsolution=optimizeCbModel(modeloATP, 'max', 'one');
Balance_energetico=ATPsolution.f

end
