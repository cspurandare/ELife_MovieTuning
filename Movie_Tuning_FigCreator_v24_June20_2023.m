%% Creating the figures for Movie tuning
clc
addpath(genpath('C:\Users\chinm\Desktop\MovieTuning\ForGithub'));
load('Allen_VVDot_DataSet_ForMovie_v1.mat')   ;

ID1=find(contains(DumpNS.Region,'LG') );
ID2=find(strcmp(DumpNS.Region,'VISp') );
ID3=find(strcmp(DumpNS.Region,'VISam') | strcmp(DumpNS.Region,'VISpm') );
ID4=find(strcmp(DumpNS.Region,'DG') );
ID5=find(strcmp(DumpNS.Region,'CA3') );
ID6=find(strcmp(DumpNS.Region,'CA1') );
ID7=find(strcmp(DumpNS.Region,'SUB') );
ID8=find(strcmp(DumpNS.Region,'APN') );
ID9=find(strcmp(DumpNS.Region,'VISlm') | strcmp(DumpNS.Region,'VISrl')...
    | strcmp(DumpNS.Region,'VISlp') | strcmp(DumpNS.Region,'VISal'));
ID10=find(strcmp(DumpNS.Region,'SCig') | strcmp(DumpNS.Region,'SCiw')...
    | strcmp(DumpNS.Region,'SCop') | strcmp(DumpNS.Region,'SCsg') | strcmp(DumpNS.Region,'SCzo'));
IDAll={ID1;ID2;ID3;ID4;ID5;ID6;ID7;};

clearvars -except DumpNS idall IDAll D

Sq1=load('Movie_BasicShf_AllTr_AllData_Seq_v2_CircSMOForSpar.mat');
Sq2=load('Movie_BasicShf_AllTr_RunFreeData_Seq_v2_CircSMOForSpar.mat');
Sq2_dash=load('Movie_BasicShf_AllTr_ShuffledRunFreeData_Seq_v2_CircSMOForSpar.mat');
Sq3=load('Movie_BasicShf_MiddleTr_AllData_Seq_v2_CircSMOForSpar.mat');
Sq5_hr=load('Movie_FRStore_Upsample_10x_seq_v7_PromScaler1p1.mat');
Sq5_hr_middleTrials=load('Movie_FRStore_Upsample_10x_seq_MiddleTrials_v7.mat');
Sq6=load('Movie_BasicShf_AllTr_SWRFreeData_Seq_v2_CircSMOForSpar.mat');
Sq6_dash=load('Movie_BasicShf_AllTr_ShuffledSWRFreeData_Seq_v2_CircSMOForSpar.mat');
Sq7=load('Movie_BasicShf_AllTr_LowerPupilArea_Seq_v2_CircSMOForSpar.mat');
Sq7_dash=load('Movie_BasicShf_AllTr_UpperPupilArea_Seq_v2_CircSMOForSpar.mat');
Sq8=load('Movie_BasicShf_AllTr_LowerThetaPower_Seq_v2_CircSMOForSpar.mat');
Sq8_dash=load('Movie_BasicShf_AllTr_UpperThetaPower_Seq_v2_CircSMOForSpar.mat');

Sc1=load('Movie_BasicShf_AllTr_AllData_Sca_v2_CircSMOForSpar.mat');
Sc5_hr=load('Movie_FRStore_Upsample_10x_sca_v7_PromScaler1p1.mat');

MV=load('Allen_natural_movies_one_and_shuffled_v2.mat');
load('Allen_ElectrodesTable_Aug2021_v2.mat');
load('Allen_SubregionStructure_Aug2021.mat')

CLR=[1 0 1 ;
    1 0 0;
    1 0.7 0;
    0 0.5 0;
    0 1 1;
    0 0 1;
    0.5882    0.4784    0.7098];
delTCurrent=0.0334;
FRCutoff_lower=0.5;
SizeBins=900;
Active=nanmean(Sq1.FR_shf_mean,2)/delTCurrent>FRCutoff_lower;
Active2=find(Active);
Active_sca=nanmean(Sc1.FR_shf_mean,2)/delTCurrent>FRCutoff_lower;
Active2_sca=find(Active);
Broad=abs(DumpNS.DUR)>0.45 & abs(DumpNS.DUR)<1.5;
Broad2=find(Broad);
Narrow=abs(DumpNS.DUR)<0.4;
Narrow2=find(Narrow);
TunedCells=find(Sq1.Zspa>2);
SpCnt2FRScaler=Sq5_hr.PRM.trHigher*Sq5_hr.PRM.delTCurrent/Sq5_hr.PRM.Upsample;
co=@colors;
c1=co('ucla gold');
c2=co('blue-violet');
c3=co('orange-red');
c4=co('ash grey');
c5=co('cyan');
ColorSeq=[0.5 0.1 0.1];
ColorSeq2=[1 0.5 0.5];
ColorSca=[0.5 0.5 1];
CLRPeak=[0 0.6 0.2];
CLRTrough=[0.6 0 0.2];
TEXTsmall={'LGN','V1','AM-PM','DG','CA3','CA1','SUB'};
A2=[27293 1416 32808 14808 26760 26735 28273;27555 29842 33374 49294 15075 14411 38138;];A2=A2';

if ~exist( 'Coh1','var')  && ~exist( 'Coh2','var')
    for i=1:900
        if i>1
            Coh1(i)=corrcoefwithNaNDiag(MV.Movie(:,:,i-1),MV.Movie(:,:,i));
            Coh2(i)=corrcoefwithNaNDiag(MV.Movie_Shf(:,:,i-1),MV.Movie_Shf(:,:,i));
        else
            Coh1(i)=corrcoefwithNaNDiag(MV.Movie(:,:,900),MV.Movie(:,:,i));
            Coh2(i)=corrcoefwithNaNDiag(MV.Movie_Shf(:,:,900),MV.Movie_Shf(:,:,i));
        end
    end
end
%% You might need to change the path for the data files you download from the Allen Brain website
% Also I kept the files in seperate folders based on the genotype of the
% mice, so the naming reflects that. If you decide on other storage
% choices, please change the directories here, accordingly
for i=1:57
    D(i).folder=strrep(D(i).folder,'L:\Matlab Codes\ProjectWise\AllenNeuropixelStuff\','**Your Directory here');
end
%% F1 Example tuned cells and CDF of tuning
%Punctate responses
%make sure A2 is defined in the first block
% A2=[1906 1416 204 12798 26760 24100 33729;27555 3215 33374 12802 27043 14411 38120];A2=A2';
%close all
figure(1);set(gcf,'Position',[748.2000  973.8000  867.2000  753.6000]);
figure(3);set(gcf,'Position',[911.4000  988.2000  867.2000  753.6000]);
UseAllData=1;
SMOVariable=[1 1 1 2 1 2 2;1 1 1 2 2 2 2;];
for j=1:size(A2,1)
    if j<4;figure(1);Subtract=1;YLIMS=[0 100];else;figure(3);Subtract=4;YLIMS=[0 50];end
    for i=1:size(A2,2)
        idcurrent=A2(j,i);
        ses=DumpNS.ExpID(idcurrent);
        FileName=[D(ses).folder '\' D(ses).name];
        UnitID=DumpNS.UnitID(idcurrent);
        FR_actC=Sq1.FR_act(idcurrent,:)/delTCurrent;
        FR_actC=jmm_smooth_1d_cor_circ(FR_actC,Sq1.PRM.SMO*SMOVariable(i,j),1);
        %FR_actC=jmm_smooth_1d_cor_circ(FR_actC,Sq1.PRM.SMO,1);
        s1=subplot(4,9,( (j-Subtract)*9 + 4*(i-1) )  + (1:4));
        set(gca,'position',get(gca,'position')+[0.00 0 -0.02 0]);
        LocIN=[s1.Position];
        if i==1 && (j-Subtract)==0
            LabelBool=1;
            annotation('textbox',[LocIN(1)-LocIN(3)*0.20 LocIN(2)+LocIN(4)*1.11 LocIN(3)/10 LocIN(4)/10],...
                'string','a','EdgeColor','none','FontSize',18,'FontWeight','bold')
        else
            LabelBool=0;
        end
        l1= Movie_PlotScatter_PSTH_psychedelic_v4_noStd(UnitID,FileName,Sq1.PRM.TCutoff,Sq1.PRM.RVCut,Sq1.PRM.FrameBins,...
            FR_actC,LocIN,LabelBool,CLR(j,:),UseAllData);
        set(gcf,'CurrentAxes',l1)
        if i==1;title([DumpNS.Region{idcurrent} '- Z_{sparsity}=' num2str(round(Sq1.Zspa(idcurrent),1))],'color',CLR(j,:));
        else;title(['Z_{sparsity}=' num2str(round(Sq1.Zspa(idcurrent),1))],'color',CLR(j,:));end
        %title([num2str(idcurrent) '-' DumpNS.Region(idcurrent)]);
        %title([num2str(idcurrent) '-' DumpNS.Region{idcurrent} '-N_{peaks}=' num2str(size(Sq5_hr.PkInfo{idcurrent},2))]);
    end
    
    p=subplot(4,9,9*(j-Subtract +1));
    set(gca,'position',get(gca,'position')+[0.00 0 0.01 0]);
    PYR=mintersect(Broad2,IDAll{j},find(~isnan(Sq1.Zspa)),Active2 );
    temp=[FracAbove2(Sq1.Zspa(PYR))  FracAbove2(Sq1.Zspa_shf(PYR)) ];
    b1=bar(1,temp(1)*100);hold on;b2=bar(2,temp(2)*100);xticks([1 2]);xticklabels({'Actual';'Chance'});
    b1.FaceColor=CLR(j,:);b1.EdgeColor=CLR(j,:);b2.FaceColor=[0.7 0.7 0.7];b2.EdgeColor=[0.7 0.7 0.7];
    if (j-Subtract)==0;ylabel('Movie tuning (%)');end
    box off; axis tight;xtickangle(25);ylim(YLIMS)
end

figure(1)
s1=subplot(4,2,7);
Zbins=linspace(-2.5,8,11);
clearvars store
for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq1.Zspa)),Active2 );store(i,1:3)=[nanlength(PYR) nansum(Sq1.Zspa(PYR)>2) FracAbove2(Sq1.Zspa(PYR))] ;
    ALL=(IDAll{i} );
    [a,b]=ecdf(Sq1.Zspa(PYR),'function','survivor');
    plot(b,a,'color',CLR(i,:),'linewidth',2)
    hold on
end
PYR=mintersect(Broad2,cell2mat(IDAll),find(~isnan(Sq1.Zspa)),Active2 );
[a,b]=ecdf(Sq1.Zspa_shf(PYR),'function','survivor');
plot(b,a,'color','k','linewidth',0.5,'linestyle','--')
hold on
xlim([-2 10]);xlabel('Movie tuning (z)')
line([2 2],ylim,'color',[0 0 0]);
legend([TEXTsmall 'Chance level']);legend boxoff;L=legend;L.ItemTokenSize(1)=8;
ylabel('Probability')
box off
grid on
LocIN=[s1.Position];
annotation('textbox',[LocIN(1)-LocIN(3)*0.20 LocIN(2)+LocIN(4)*1.11 LocIN(3)/10 LocIN(4)/10],...
    'string','d','EdgeColor','none','FontSize',18,'FontWeight','bold')

s1=subplot(4,2,8);
hold on;clearvars b
for i=1:7;b(i)=bar(i,100*store(i,3));b(i).FaceColor=CLR(i,:);b(i).EdgeColor=CLR(i,:);end
xticks([1:7]);xticklabels(TEXTsmall);axis tight;ylabel('Movie tuning (%)')
LocIN=[s1.Position];
annotation('textbox',[LocIN(1)-LocIN(3)*0.20 LocIN(2)+LocIN(4)*1.11 LocIN(3)/10 LocIN(4)/10],...
    'string','e','EdgeColor','none','FontSize',18,'FontWeight','bold')
xtickangle(25)
hold on
line(xlim,[2.25 2.25],'color',[0 0 0],'linestyle','--')
%% ED Plot the same neurons for stationary only data, and compare fraction selective in all, runfree and equi-subsample
figure('position',[510.6000  -22.2000  989.6000  775.2000]);
SPLoc=[1:7];
UseAllData=0;
A2Current=A2;
A2Current(4,1)=27048;
for j=1:size(A2Current,1)
    for i=1:1%size(A2,2)
        idcurrent=A2Current(j,i);
        ses=DumpNS.ExpID(idcurrent);
        FileName=[D(ses).folder '\' D(ses).name];
        UnitID=DumpNS.UnitID(idcurrent);
        FR_actC=Sq2.FR_act(idcurrent,:)/delTCurrent;
        FR_actC=jmm_smooth_1d_cor_circ(FR_actC,Sq1.PRM.SMO*SMOVariable(i,j),1);
        %FR_actC=jmm_smooth_1d_cor_circ(FR_actC,Sq1.PRM.SMO,1);
        s1=subplot(3,3,SPLoc(j));
        LocIN=[s1.Position];if i==1; LocIN= LocIN+[0 0.0 0 0];end
        if  j==1
            LabelBool=1;
            annotation('textbox',[LocIN(1)-LocIN(3)*0.20 LocIN(2)+LocIN(4)*1.11 LocIN(3)/10 LocIN(4)/10],...
                'string','A','EdgeColor','none','FontSize',18,'FontWeight','bold')
        else
            LabelBool=0;
        end
        l1= Movie_PlotScatter_PSTH_psychedelic_v4_noStd(UnitID,FileName,Sq2.PRM.TCutoff,Sq2.PRM.RVCut,Sq2.PRM.FrameBins,...
            FR_actC,LocIN,LabelBool,CLR(j,:),UseAllData);
        set(gcf,'CurrentAxes',l1)
        title(TEXTsmall{j},'color',CLR(j,:));
        %title([num2str(idcurrent) '-' DumpNS.Region(idcurrent)]);
    end
end

Zbins=linspace(-2.5,8,11);
clearvars store
s1=subplot(3,3,8);
for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq2.Zspa)),Active2 );store(i,1:3)=[nanlength(PYR) nansum(Sq2.Zspa(PYR)>2) FracAbove2(Sq2.Zspa(PYR))] ;
    ALL=(IDAll{i} );
    [a,b]=ecdf(Sq2.Zspa(PYR),'function','survivor');
    plot(b,a,'color',CLR(i,:),'linewidth',2)
    hold on
end
PYR=mintersect(Broad2,cell2mat(IDAll),find(~isnan(Sq2.Zspa)),Active2 );
[a,b]=ecdf(Sq2.Zspa_shf(PYR),'function','survivor');
plot(b,a,'color','k','linewidth',0.5,'linestyle','--')
hold on
xlim([-2 10]);xlabel('Movie tuning (z)','FontSize',8)
line([2 2],ylim,'color',[0 0 0]);
legend(TEXTsmall);legend boxoff;L=legend;L.ItemTokenSize(1)=8;
ylabel('Probability')
box off
grid on
LocIN=[s1.Position];
annotation('textbox',[LocIN(1)-LocIN(3)*0.20 LocIN(2)+LocIN(4)*1.11 LocIN(3)/10 LocIN(4)/10],...
    'string','B','EdgeColor','none','FontSize',18,'FontWeight','bold')

clearvars storeAll storeSta storeEq ks_pvalue
for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq1.Zspa)),Active2 );
    storeAll(i,1:3)=[nanlength(PYR) nansum(Sq1.Zspa(PYR)>2) FracAbove2(Sq1.Zspa(PYR))] ;
    PYR2d=mintersect(Broad2,IDAll{i},find(~isnan(Sq2.Zspa)),Active2 );
    storeSta(i,1:3)=[nanlength(PYR2d) nansum(Sq2.Zspa(PYR2d)>2) FracAbove2(Sq2.Zspa(PYR2d))] ;
    PYR2d=mintersect(Broad2,IDAll{i},find(~isnan(Sq2_dash.Zspa)),Active2 );
    storeEq(i,1:3)=[nanlength(PYR2d) nansum(Sq2_dash.Zspa(PYR2d)>2) FracAbove2(Sq2_dash.Zspa(PYR2d))] ;
    [~,ks_pvalue(i)]=kstest2(Sq2.Zspa(PYR),Sq2_dash.Zspa(PYR2d));
end
s1=subplot(3,3,9);
b=bar([storeAll(:,end) storeSta(:,end) storeEq(:,end)]*100);
box off;
for i=1:3
    b(i).EdgeColor=b(i).FaceColor;
end
box off
ylabel('Movie tuning (%)')
xticks([1:7])
xticklabels(TEXTsmall)
xtickangle(20)
line(xlim,[2.25 2.25],'color',[0 0 0],'linestyle','--')
legend('All data','Stationary data','Equivalent subsample','Chance level')
L=legend;legend boxoff;L.ItemTokenSize(1)=8;
LocIN=[s1.Position];
annotation('textbox',[LocIN(1)-LocIN(3)*0.20 LocIN(2)+LocIN(4)*1.11 LocIN(3)/10 LocIN(4)/10],...
    'string','C','EdgeColor','none','FontSize',18,'FontWeight','bold')
%% F2 Peak and trough stuff - This based on the latest upsampled tuning curves, and I have found the
% peaks already and saved them in the Sq5_HR above
%close all
BoolExtractForFig3=0;
clearvars p_chisq SummaryPeakProp DumpPeakProp DumpPeakProp2 SummaryPeakProp XCorPeakHist CorrCf
pd=makedist('Uniform',1,900);
for i=1:7
    clearvars TC sTC tc2 stc2
    idc=mintersect(Broad2,IDAll{i},TunedCells,Active2);
    TC=Sq5_hr.SpCnt_actual(idc,:);
    sTC=Sq5_hr.SpCnt_shuffle(idc,:);
    for sm=1:10
        tc2(sm,:,:)=colsmoothChinmay_Circ(TC,Sq5_hr.PRM.SmoVec(sm)*Sq5_hr.PRM.Upsample);
        stc2(sm,:,:)=colsmoothChinmay_Circ(sTC,Sq5_hr.PRM.SmoVec(sm)*Sq5_hr.PRM.Upsample);
    end
    PKI=Sq5_hr.PkInfo(idc);
    PK_SpCnt=Sq5_hr.SpCnt_Vect_Trialwise(idc);
    binCenters=center(Sq5_hr.PRM.FrameBins);
    
    clearvars pkCount pkWidths pkScores pkWidthRatio pkMedianWidth pkTotWidth pkMedianArea pkTotArea pkAreas pkLocs pk_spikingCOV0 pk_meanRate CellIds SesIds
    pkWidths=[];
    pkScores=[];
    pkAreas=[];
    pkLocs=[];
    pk_spikingCOV0=[];
    pk_meanRate=[];
    CellIds=[];
    SesIds=[];
    
    PkExtract=cell(size(PKI,1),1);
    for j=1:length(PKI)
        pkCount(j)=size(PKI{j},2);
        
        PkArea=nan(size(PKI{j},2),2);
        for pk=1:size(PKI{j},2)
            BinsInterest=isbetween(binCenters,PKI{j}(2,pk)-PKI{j}(3,pk),PKI{j}(2,pk)+PKI{j}(4,pk) );
            PkArea(pk,:)=[nansum(TC(j,BinsInterest)) ; ...
                (nansum(sTC(j,:)) * nansum(BinsInterest)/nanlength(BinsInterest) )] ;
            if BoolExtractForFig3==1
                PkExtract{j}(pk,1:size(TC,2))=nan;
                tc_curr=tc2(PKI{j}(7,pk)==Sq5_hr.PRM.SmoVec,j,:);
                stc_curr=stc2(PKI{j}(7,pk)==Sq5_hr.PRM.SmoVec,j,:);
                PkExtract{j}(pk,BinsInterest)=tc_curr(BinsInterest)./nanmean(stc_curr);
            end
        end
        
        if size(PKI{j},2)>0
            
            if size(PKI{j},2)>1
                pkWidthRatio(j)=nanmax(Sq1.PRM.delTCurrent*(PKI{j}(3,:) + PKI{j}(4,:)))/...
                    nanmin(Sq1.PRM.delTCurrent*(PKI{j}(3,:) + PKI{j}(4,:)));
            else
                pkWidthRatio(j)=nan;
            end
            pkTotWidth(j)=Sq1.PRM.delTCurrent*(nansum(PKI{j}(3,:) + PKI{j}(4,:)));
            pkTotArea(j)=nansum(PkArea(:,1))/ nansum(PkArea(:,2));
            
            pkMedianArea(j)=nanmedian(PkArea(:,1)./PkArea(:,2));
            pkMedianWidth(j)=Sq1.PRM.delTCurrent*(nanmedian(PKI{j}(3,:) + PKI{j}(4,:)));
            
            pkWidths=[pkWidths Sq1.PRM.delTCurrent*(PKI{j}(3,:) + PKI{j}(4,:))];
            pkScores=[pkScores PKI{j}(6,:)];
            pkAreas=[pkAreas (PkArea(:,1)./PkArea(:,2))'];
            CellIds=[CellIds ones(1,size(PKI{j},2))*idc(j)];
            SesIds=[SesIds DumpNS.ExpID(idc(j))*ones(1,size(PKI{j},2));];
            
            pkLocs=[pkLocs PKI{j}(2,:) ];
            %pk_spikingCOV0=[pk_spikingCOV0 nanstd(PK_SpCnt{j}')./nanmean(PK_SpCnt{j}');];
            %pk_meanRate=[pk_meanRate nanmean(PK_SpCnt{j}')./ (Sq1.PRM.delTCurrent*(PKI{j}(3,:) + PKI{j}(4,:)))];
        else
            pkMedianWidth(j)=nan;
            pkWidthRatio(j)=nan;
            pkTotWidth(j)=nan;
            pkMedianArea(j)=nan;
            pkTotArea(j)=nan;
        end
    end
    DumpPeakProp{i,1}=[pkCount];DumpPeakProp{i,2}=[pkWidthRatio];DumpPeakProp{i,3}=[pkTotWidth];DumpPeakProp{i,4}=[pkTotArea];
    DumpPeakProp{i,5}=[pkWidths];DumpPeakProp{i,6}=[pkAreas]; DumpPeakProp{i,7}=[pkMedianWidth];DumpPeakProp{i,8}=[pkMedianArea];
    DumpPeakProp{i,9}=[CellIds];DumpPeakProp{i,10}=[SesIds];DumpPeakProp{i,11}=cell2mat(PkExtract);DumpPeakProp{i,12}=pkLocs;
    for summ=1:8
        SummaryPeakProp(i,summ,:)=[nanmean(DumpPeakProp{i,summ}) nanmedian(DumpPeakProp{i,summ})...
            nanstd(DumpPeakProp{i,summ}) nanstd(DumpPeakProp{i,summ})/sqrt(nanlength(DumpPeakProp{i,summ})) ...
            nanstd(DumpPeakProp{i,summ})/sqrt(nanlength(DumpPeakProp{i,summ}))*1.96];
    end
    DumpTemp{i}=[pkWidths;pkScores];
    %__________________________________________________________
    figure(1);if i==7;set(gcf,'position',[330.6000   15.4000  944.8000  750.4000]);delete(findall(gcf,'type','annotation'));end
    %Number of peaks, all widths, median widths, total frames under peaks,
    %total area under peaks
    subplot_CSP(3,3,1,1)
    hold on
    [a,b]=ecdf(pkCount);plot(b,a,'linewidth',2,'color',CLR(i,:));
    hold on
    if i==7;xlabel('# Peaks');grid on;xlim([0 20]);box off;L=legend(TEXTsmall);L.ItemTokenSize(1)=7;legend boxoff;ylabel('Probability');end
    
    subplot_CSP(3,3,4,2);
    hold on
    [a,b]=ecdf(pkWidths);plot(b,a,'linewidth',2,'color',CLR(i,:));
    hold on
    if i==7;xlabel('All peak widths (sec)');grid on;xlim([0.01 25]);xticks([0.01 0.1 1 10 25]);box off;set(gca,'xscale','log');    end
    
    subplot_CSP(4,3,10,3);
    hold on
    [a,b]=ecdf(pkTotWidth);plot(b,a,'linewidth',2,'color',CLR(i,:));nanmean(pkTotWidth)
    hold on
    if i==7;xlabel('Cumulative width of peaks (sec)');grid on;xlim([0 25]);xticks([0 10 20]);box off;end
    
    subplot_CSP(4,3,11,4);
    hold on
    [a,b]=ecdf(pkWidthRatio);plot(b,a,'linewidth',2,'color',CLR(i,:));
    hold on
    if i==7;xlabel('Peak width ratio - Broad/Narrow');grid on;xlim([1 1000]);xticks([1 10 100]);box off;set(gca,'xscale','log');end
    
    subplot_CSP(4,3,12,5);
    hold on
    [a,b]=ecdf(pkAreas);plot(b,a,'linewidth',2,'color',CLR(i,:));
    hold on
    if i==7;L=legend(TEXTsmall) ;legend boxoff;L.ItemTokenSize(1)=8;xlabel('All firing in peaks (norm.)');grid on;xlim([1 20]);xticks([1 2 4 8 16]);box off;set(gca,'xscale','log');end
    
    %___________________________Additional things
    figure(2);if i==7;set(gcf,'position',1e3*[0.1706    0.0898    1.1968    0.5456]);delete(findall(gcf,'type','annotation'));end
    
    subplot_CSP(3,3,1)
    hold on
    [a,b]=ecdf(pkMedianWidth);plot(b,a,'linewidth',2,'color',CLR(i,:));
    hold on
    if i==7;xlabel('Median peak width (sec)');grid on;xlim([0.02 25]);xticks([ 0.1 1 10 25]);box off;set(gca,'xscale','log');   end
    
    subplot_CSP(3,3,2);
    hold on
    [a,b]=ecdf(pkMedianArea);plot(b,a,'linewidth',2,'color',CLR(i,:));
    hold on
    if i==7;xlabel('Median firing in peaks (norm.)');grid on;xlim([1 8]);xticks([1 2 4 8]);box off;set(gca,'xscale','log');end
    
    subplot_CSP(3,3,3);
    hold on
    [a,b]=ecdf(pkTotArea);plot(b,a,'linewidth',2,'color',CLR(i,:));
    hold on
    if i==7;xlabel('Cumulative firing in peaks (norm.)');grid on;xlim([1 8]);xticks([1 2 4 8]);box off;set(gca,'xscale','log');end
    
    clearvars RatioAct RatioShf
    SesAvail=unique(SesIds);
    for ses=SesAvail
        CurrentID=SesIds==ses;
        CurrentWidths=pkWidths(CurrentID);
        CurrentCellID=CellIds(CurrentID);
        for n=1:1
            CellIDShf=CurrentCellID(randperm(length(CurrentCellID)));
            CellOptions=unique(CurrentCellID);
            for un=1:length(CellOptions)
                wdths_act=CurrentWidths(CurrentCellID==CellOptions(un));
                wdths_shf=CurrentWidths(CellIDShf==CellOptions(un));
                if length(wdths_act)>=2
                    RatioAct(idc==CellOptions(un))=max(wdths_act)/min(wdths_act);
                    RatioShf(idc==CellOptions(un))=max(wdths_shf)/min(wdths_shf);
                else
                    RatioAct(idc==CellOptions(un))=nan;
                    RatioShf(idc==CellOptions(un))=nan;
                end
            end
        end
    end
    if i==1
        subplot_CSP(3,7,i+7,4);
    else
        subplot(3,7,i+7);
    end
    [a,b]=ecdf(RatioAct);plot(b,a,'linewidth',2,'color',CLR(i,:))
    hold on
    [a,b]=ecdf(RatioShf);plot(b,a,'linewidth',2,'color',[0.7 0.7 0.7])
    set(gca,'xscale','log');
    [~,RatioKS_pval(i)]=kstest2(log10(RatioAct),log10(RatioShf));
    DumpRatios{i}=[RatioAct;RatioShf];
    title(TEXTsmall{i},'color',CLR(i,:));box off
    xlim([1 1e3]);xticks([1 10 100 1000])
    if i==1
        L=legend('Actual','Shuffle grouping');xlabel('Widest/Narrowest - Peak width ratio')
        legend boxoff;L.ItemTokenSize(1)=7;ylabel('Probability');xticks([1 10 100]);
    end
    %_________________________Width vs peak score stuff now
    if i==1
        subplot_CSP(3,7,i+14,5);
    else
        subplot(3,7,i+14);
    end
    binWidths=linspace(-2,1.4,30);
    binScore=linspace(0.04,1,30);
    DensityCount=hist3(log10(DumpTemp{i}'),{binWidths binScore});
    DensityCount=100*smooth2a(DensityCount,1)/nansum(DensityCount(:));
    DensityCount(DensityCount==0)=nan;
    DensityCount=log10(DensityCount);
    imagescwithnan(center(binWidths), center(binScore), DensityCount',jet,[1 1 1]);
    caxis([-3.5 0])
    if i==1
        cbh = colorbar ;
        cbh.Ticks=[-3 -2 -1 0 1];
        %cbh.TickLabels = num2cell(roundsd(10.^cbh.Ticks,2));
        ylabel(cbh,'Peaks (%)','FontSize',10);
    end
    fg=gca;fg.XTickLabel = num2cell(roundsd(10.^fg.XTick,2));fg.YTickLabel = num2cell(roundsd(10.^fg.YTick,2));
    %cla;plot(DumpTemp{i}(1,:),DumpTemp{i}(2,:),'.','color',CLR(i,:));set(gca,'xscale','log');set(gca,'yscale','log');
    title(TEXTsmall{i},'color',CLR(i,:));axis square;
    if i==1;xlabel('Peak width (sec)');ylabel('Peak prominence (norm.)');title([TEXTsmall{i} ],'color',CLR(i,:));end
    
   
    if i<4
        idc_sca=mintersect(idc,Broad2,IDAll{i},find(Sc1.Zspa>2),Active2_sca);
        PKI=Sc5_hr.PkInfo(idc_sca);
        TC=Sc5_hr.SpCnt_actual(idc_sca,:);
        sTC=Sc5_hr.SpCnt_shuffle(idc_sca,:);
        
        clearvars pkCount_sca pkWidths_sca pkWidthRatio_sca pkMedianWidth_sca pkTotWidth_sca pkMedianArea_sca pkTotArea_sca pkAreas pkLocs_sca CellIds_sca
        
        pkWidths_sca=[];
        pkAreas_sca=[];
        pkLocs_sca=[];
        CellIds_sca=[];
        for j=1:length(PKI)
            pkCount_sca(j)=size(PKI{j},2);
            
            PkArea=nan(size(PKI{j},2),2);
            for pk=1:size(PKI{j},2)
                BinsInterest=isbetween(binCenters,PKI{j}(2,pk)-PKI{j}(3,pk),PKI{j}(2,pk)+PKI{j}(4,pk) );
                PkArea(pk,:)=[nansum(TC(j,BinsInterest)) ; ...
                    (nansum(sTC(j,:)) * nansum(BinsInterest)/nanlength(BinsInterest) )] ;
            end
            
            if size(PKI{j},2)>0
                
                if size(PKI{j},2)>1
                    pkWidthRatio_sca(j)=nanmax(Sq1.PRM.delTCurrent*(PKI{j}(3,:) + PKI{j}(4,:)))/...
                        nanmin(Sq1.PRM.delTCurrent*(PKI{j}(3,:) + PKI{j}(4,:)));
                else
                    pkWidthRatio_sca(j)=nan;
                end
                pkTotWidth_sca(j)=Sq1.PRM.delTCurrent*(nansum(PKI{j}(3,:) + PKI{j}(4,:)));
                pkTotArea_sca(j)=nansum(PkArea(:,1))/ nansum(PkArea(:,2));
                
                pkMedianArea_sca(j)=nanmedian(PkArea(:,1)./PkArea(:,2));
                pkMedianWidth_sca(j)=Sq1.PRM.delTCurrent*(nanmedian(PKI{j}(3,:) + PKI{j}(4,:)));
                
                pkWidths_sca=[pkWidths_sca Sq1.PRM.delTCurrent*(PKI{j}(3,:) + PKI{j}(4,:))];
                pkAreas_sca=[pkAreas_sca (PkArea(:,1)./PkArea(:,2))'];
                CellIds_sca=[CellIds_sca ones(1,size(PKI{j},2))*idc(j)];
                
                pkLocs_sca=[pkLocs_sca PKI{j}(2,:) ];
                
            else
                pkMedianWidth_sca(j)=nan;
                pkWidthRatio_sca(j)=nan;
                pkTotWidth_sca(j)=nan;
                pkMedianArea_sca(j)=nan;
                pkTotArea_sca(j)=nan;
            end
        end
        DumpPeakProp2{i,1}=[pkCount_sca];DumpPeakProp2{i,2}=[pkWidthRatio_sca];DumpPeakProp2{i,3}=[pkTotWidth_sca];DumpPeakProp2{i,4}=[pkTotArea_sca];
        DumpPeakProp2{i,5}=[pkWidths_sca];DumpPeakProp2{i,6}=[pkAreas_sca]; DumpPeakProp2{i,7}=[pkMedianWidth_sca];DumpPeakProp2{i,8}=[pkMedianArea_sca];
        DumpPeakProp2{i,9}=[CellIds_sca];
        for summ=1:8
            SummaryPeakProp2(i,summ,:)=[nanmean(DumpPeakProp2{i,summ}) nanmedian(DumpPeakProp2{i,summ})...
                nanstd(DumpPeakProp2{i,summ}) nanstd(DumpPeakProp2{i,summ})/sqrt(nanlength(DumpPeakProp2{i,summ})) ...
                nanstd(DumpPeakProp2{i,summ})/sqrt(nanlength(DumpPeakProp2{i,summ}))*1.96];
        end
        
        figure(4);if i==3;set(gcf,'position',1e3*[0.2186    0.1330    1.0000    0.5000]);delete(findall(gcf,'type','annotation'));end
        %For only the visual areas, repeat peak quantification
        %Number of peaks, all widths, median widths, total frames under peaks,
        %total area under peaks
        subplot_CSP(2,3,1);
        hold on
        [a,b]=ecdf(pkCount_sca);plot(b,a,'linewidth',2,'color',CLR(i,:));
        hold on
        if i==3;xlabel('# Peaks');grid on;xlim([0 85]);box off;L=legend(TEXTsmall(1:3));L.ItemTokenSize(1)=7;legend boxoff;end
        
        subplot_CSP(2,3,2);
        hold on
        [a,b]=ecdf(pkMedianWidth_sca(~isinf(pkMedianWidth_sca)));plot(b,a,'linewidth',2,'color',CLR(i,:));
        hold on
        if i==3;xlabel('Median peak width (sec)');grid on;xlim([0.02 10]);box off;set(gca,'xscale','log');   end
        
        subplot_CSP(2,3,3);
        hold on
        [a,b]=ecdf(pkWidths_sca(~isinf(pkWidths_sca)));plot(b,a,'linewidth',2,'color',CLR(i,:));
        hold on
        if i==3;xlabel('Peak widths (sec)');grid on;xlim([0.015 10]);box off;set(gca,'xscale','log');    end
        
        subplot_CSP(2,3,4);
        hold on
        [a,b]=ecdf(pkTotWidth_sca(~isinf(pkTotWidth_sca)));plot(b,a,'linewidth',2,'color',CLR(i,:));
        hold on
        if i==3;xlabel('Cumulative width of peaks (sec)');grid on;xlim([0 10]);box off;end
        
        subplot_CSP(2,3,5);
        hold on
        [a,b]=ecdf(pkWidthRatio_sca(~isinf(pkMedianArea_sca)));plot(b,a,'linewidth',2,'color',CLR(i,:));
        hold on
        if i==3;xlabel('Peak width ratio - Broad/Narrow');grid on;xlim([1 1000]);box off;set(gca,'xscale','log');end
        
        subplot_CSP(2,3,6);
        hold on
        [a,b]=ecdf(pkTotArea_sca(~isinf(pkTotArea_sca)));plot(b,a,'linewidth',2,'color',CLR(i,:));
        hold on
        if i==3;xlabel('Total firing in peaks (norm.)');grid on;xlim([1 20]);box off;set(gca,'xscale','log');end
        
        figure(5);% compare widths and number of peaks
        set(gcf,'position',1e3*[ 0.2882    0.0498    1.0000    0.7328])
        clearvars matchBack
        for index=1:length(idc_sca)
            matchBack(index)=find(idc==idc_sca(index));
        end
        if i==1;subplot_CSP(4,6,1+6*(i-1),1);else;subplot(4,6,1+6*(i-1));end
        plot(pkCount(matchBack)+0.5-rand(size(pkCount(matchBack))),pkCount_sca+0.5-rand(size(pkCount(matchBack))),'.','color',CLR(i,:),'markersize',2)
        xlim([0 70]);ylim([0 70]);box off;axis square
        if i==1;xlabel('# Movie-fields - continuous');ylabel('Scrambled');end
        plot45Line([0 0 0 ])
        title(TEXTsmall{i},'color',CLR(i,:));
        
        if i==1;subplot_CSP(4,6,2+6*(i-1),2);else;subplot(4,6,2+6*(i-1));end
        plot(pkMedianWidth(matchBack),pkMedianWidth_sca,'.','color',CLR(i,:),'markersize',2)
        xlim([0.1 30]);ylim([0.1 30]);xticks([0.1 1 10]);yticks([0.1 1 10])
        plot45Line([0 0 0 ])
        set(gca,'xscale','log');set(gca,'yscale','log');
        if i==1;xlabel('Median duration - continuous (sec)');ylabel('Scrambled (sec)');end
        box off;axis square
        
        if i==1;subplot_CSP(4,6,3+6*(i-1),3);else;subplot(4,6,3+6*(i-1));end
        plot(pkTotWidth(matchBack),pkTotWidth_sca,'.','color',CLR(i,:),'markersize',2)
        xlim([0.1 30]);ylim([0.1 30]);box off;axis square
        plot45Line([0 0 0 ])
        if i==1;xlabel('Cumulative duration - continuous (sec)');ylabel('Scrambled (sec)');end
        box off;axis square
        %__________________Histograms now__________________________________
        PkCountBins=[0:2:70];
        if i==1;subplot_CSP(4,6,4+6*(i-1),4);else;subplot(4,6,4+6*(i-1));end
        [a,b]=histperc(pkCount(matchBack)+0.5-rand(size(pkCount(matchBack))),PkCountBins);
        stairs(b,jmm_smooth_1d_cor(a,2,1),'linewidth',2,'color',ColorSeq);
        hold on
        [a,b]=histperc(pkCount_sca+0.5-rand(size(pkCount(matchBack))),PkCountBins);
        stairs(b,jmm_smooth_1d_cor(a,2,1),'linewidth',2,'color',ColorSca);
        box off;axis tight;xticks([0 30 60]);grid on;ylim([0 15]);yticks([0 15])
        if i==1;xlabel('# Fields');ylabel('Cells (%)');L=legend('Continuous','Scrambled');legend boxoff;L.ItemTokenSize(1)=8;end
        
        LogBins=logspace(-2.1,1.3,40);
        LinBins=linspace(0,20,40);
        if i==1;subplot_CSP(4,6,5+6*(i-1),5);else;subplot(4,6,5+6*(i-1));end
        [a,b]=histperc(pkMedianWidth(matchBack),LogBins);
        stairs(b,jmm_smooth_1d_cor(a,2,1),'linewidth',2,'color',ColorSeq);
        hold on
        [a,b]=histperc(pkMedianWidth_sca,LogBins);
        stairs(b,jmm_smooth_1d_cor(a,2,1),'linewidth',2,'color',ColorSca);
        box off;axis tight;ylim([0 16]);xticks([0.01 0.1 1 10]);grid on
        if i==1;xlabel('Median field duration (sec)');end
        set(gca,'xscale','log');
        
        if i==1;subplot_CSP(4,6,6+6*(i-1),6);else;subplot(4,6,6+6*(i-1));end
        [a,b]=histperc(pkTotWidth(matchBack),LinBins);
        stairs(b,jmm_smooth_1d_cor(a,2,1),'linewidth',2,'color',ColorSeq);
        hold on
        [a,b]=histperc(pkTotWidth_sca,LinBins);
        stairs(b,jmm_smooth_1d_cor(a,2,1),'linewidth',2,'color',ColorSca);
        box off;axis tight;ylim([0 16]);xticks([ 5 10 15 20]);grid on
        if i==1;xlabel('Cumulative duration (sec)');end
        %---------------------------Ratio now-------------------
        subplot_CSP(4,3,10,7);
        Ratio=( pkCount(matchBack)+0.5-rand(size(pkCount(matchBack))) )./ ( pkCount_sca+0.5-rand(size(pkCount(matchBack))) );
        [a,b]=ecdf(Ratio);plot(b,a,'color',CLR(i,:),'linewidth',1.5);
        hold on
        if i==1;ylabel('Probability');xlabel('# Fields ratio - Continuous / Scrambled');end
        set(gca,'xscale','log');xlim([0.1 100]);box off;grid on
        DumpSeqScaRatio{i,1}=Ratio;
        
        subplot_CSP(4,3,11,8);
        Ratio=( pkMedianWidth(matchBack) )./ ( pkMedianWidth_sca );
        [a,b]=ecdf(Ratio);plot(b,a,'color',CLR(i,:),'linewidth',1.5);
        hold on
        if i==1;xlabel('Median duration ratio ');end
        set(gca,'xscale','log');xlim([0.1 100]);box off;grid on
        DumpSeqScaRatio{i,2}=Ratio;
        
        subplot_CSP(4,3,12,9);
        Ratio=( pkTotWidth(matchBack) )./ ( pkTotWidth_sca );
        [a,b]=ecdf(Ratio);plot(b,a,'color',CLR(i,:),'linewidth',1.5);
        hold on
        if i==1;xlabel('Cumulative duration ratio');end
        xlim([0.1 20]);box off;grid on
        DumpSeqScaRatio{i,3}=Ratio;
        
        DumpPeakProp3{i,1}=pkCount(matchBack);
        DumpPeakProp3{i,2}=pkTotWidth(matchBack);
        DumpPeakProp3{i,3}=pkMedianWidth(matchBack);
        
        for summ=1:3
            SummaryPeakProp3(i,summ,:)=[nanmean(DumpPeakProp3{i,summ}) nanmedian(DumpPeakProp3{i,summ})...
                nanstd(DumpPeakProp3{i,summ}) nanstd(DumpPeakProp3{i,summ})/sqrt(nanlength(DumpPeakProp3{i,summ})) ...
                nanstd(DumpPeakProp3{i,summ})/sqrt(nanlength(DumpPeakProp3{i,summ}))*1.96];
        end
    end
    
end
if BoolExtractForFig3==1
    save('Movie_PeakRelatedData_DumpPeakProps_June24_2022.mat','DumpPeakProp','DumpPeakProp2','DumpPeakProp3','SummaryPeakProp','SummaryPeakProp2','-v7.3');
end

figure(1)
CellIds=[3209   ;3209   ;26754     ;26754  ;27286  ;32803;
    43167;  54050;   38467;   28437;     48906;   43285;];
Limits=[153 156 ;650 750; 194 197.5;440 540;550 850;  0 6;
    67 73;359 365; 182 192; 462 468; 708 710.5; 136 144;];
PlotLoc=[3 4 5 6 9 10 11 12 15:18];
SMOCurrent=[2 128 4 128 128 4 4*ones(1,6)];
LabelBool=[1 0 0 0 0 0 zeros(1,6)];
ScatterSize=[4 4 6 2 2 4 4*ones(1,6)];
RegionCurrent=[2 2 6 6 1 3 4 5 6 7 7 4];

TTL1=0.05/Sq1.PRM.delTCurrent;
TTL2=1/Sq1.PRM.delTCurrent;
TimeTickLength=[TTL1 TTL2 TTL1 TTL2 TTL2 TTL1*ones(1,7)];
for i=1:12
    idcurrent=CellIds(i);
    ses=DumpNS.ExpID(idcurrent);
    FileName=[D(ses).folder '\' D(ses).name];
    UnitID=DumpNS.UnitID(idcurrent);
    FR_actC=jmm_smooth_1d_cor_circ(Sq5_hr.SpCnt_actual(idcurrent,:)/SpCnt2FRScaler,SMOCurrent(i),1);
    PksAccepted=Sq5_hr.PkInfo{idcurrent};
    FR_at_peaks=FR_actC(nearestpoint(PksAccepted(2,:),center(Sq5_hr.PRM.FrameBins)));
    switch i
        case 1
            subplot_CSP(4,6,PlotLoc(i),6);
        case 3
            subplot_CSP(4,6,PlotLoc(i),7);
        case 5
            subplot_CSP(4,6,PlotLoc(i),8);
        otherwise
            subplot(4,6,PlotLoc(i));
    end
    
    Movie_PlotScatter_psychedelic_v5_hr_Overlapping(UnitID,FileName,Sq1.PRM.TCutoff,Sq1.PRM.RVCut,LabelBool(i),CLR(RegionCurrent(i),:),ScatterSize(i),FR_actC,center(Sq5_hr.PRM.FrameBins));
    xlim([Limits(i,1) Limits(i,2)])
    XMaxTBU=Limits(i,2)-0.5;
    yyaxis left;
    line([XMaxTBU (XMaxTBU-TimeTickLength(i))],[63 63],'color','k','linewidth',2);
    line([XMaxTBU XMaxTBU],[60.5 63],'color','k','linewidth',2);
    line([(XMaxTBU-TimeTickLength(i)) (XMaxTBU-TimeTickLength(i))],...
        [60.5 63],'color','k','linewidth',2);
    ylim([1 63.5]);
    yyaxis right;YL=ylim;ylim([YL(1) YL(2)*63.5/60])
end

%Check if the cumulative width is similar for CA1 and V1 by chance
temp1=[DumpPeakProp{2,5} ;DumpPeakProp{2,9}];Uni1=unique(temp1(2,:));temp1=temp1(:,1:find(temp1(2,:)==Uni1(2000)));
temp2=[DumpPeakProp{6,5} ;DumpPeakProp{6,9}];Uni2=unique(temp2(2,:));temp2=temp2(:,1:find(temp2(2,:)==Uni2(2000)));

%I have peaks from 2000 tuned neurons, stored out here, so now I can
%compare cumulative width, actual and shuffled
clearvars CumWidth_V1 CumWidth_CA1
for k=1:2000
    CumWidth_V1(k,1)=nansum(temp1(1,temp1(2,:)==Uni1(k)));%How many peaks for this neuron, actual
    CumWidth_V1(k,2)=nansum(temp1(1,randperm(size(temp1,2) , sum(temp2(2,:)==Uni2(k)) ) ));
    % For this neuron #, in CA1, check how many peaks it had, and then get
    % that many widths to add up, from the V1 widths
    CumWidth_CA1(k,1)=nansum(temp2(1,temp2(2,:)==Uni2(k)));
    CumWidth_CA1(k,2)=nansum(temp2(1,randperm(size(temp2,2) , sum(temp1(2,:)==Uni1(k)) ) ));
end

%% F3 Stack plot, MSUA, and peak trough counts
%close all
clearvars AvgNonDiag
SqTBU=Sq1;
load('Movie_PeakRelatedData_DumpPeakProps_May12_2022.mat');
binCenters=center(Sq5_hr.PRM.FrameBins);
BinFrame=linspace(binCenters(1),binCenters(end),60);
pd=makedist('Uniform',1,900);
for i=1:7
    %p=subplot(8,3,1+(i-1)*3);
    figure(1)
    set(gcf,'position',1e3*[0.2418    0.1162    1.0496    0.6256]);
    %-------------------------------Peak stack-------------------
    if i==1;subplot_CSP(4,7,i,1);else;subplot(4,7,i);end
    %SnakePlotMaker(center(Sq5_hr.PRM.FrameBins),colNormalizem1p1(DumpPeakProp{i,11}));
    if i==1;xlabel('Frame #');ylabel('Movie-field #');
    end
    box off
    xticks([400 800]);yticks([0 floor(size(DumpPeakProp{i,11},1)/10)*10]);
    title(TEXTsmall{i},'color',CLR(i,:));
    %------------------------------Snake plot ----------------------
    if i==1;subplot_CSP(4,7,i+7,2);else;subplot(4,7,i+7);end
    PYR=mintersect(find(SqTBU.Zspa>2),IDAll{i},Broad2,Active2 );
    TCs=SqTBU.FR_act(PYR,:);
    TCs=colmaxNormalize(colsmoothChinmay_Circ(TCs,Sq1.PRM.SMO*5));
    SnakePlotMaker_removeNan(center(SqTBU.PRM.FrameBins),TCs);
    if i==1;xlabel('Frame #');ylabel('Cell #');
    end
    colorbar off;axis tight;ylim([0 inf]);%caxis([0 1])
    xticks([400 800]);yticks([0 floor(size(TCs(~isnan(TCs(:,1)),:),1)/10)*10])
    box off
    %title(TEXTsmall{i},'color',CLR(i,:))
    ax=gca;ax.YAxis(1).Color = CLR(i,:);ax.XAxis(1).Color = CLR(i,:);
    temp=TCs(:);
    AvgNonDiag(i)=nanmean(temp(temp~=1));
    %-------------------------------MSUA response-----------------
    if i==1;subplot_CSP(4,7,i+14,3);else;subplot(4,7,i+14);end
    TCs=SqTBU.FR_act(PYR,:);
    TCs=colmaxNormalize(colsmoothChinmay_Circ(TCs,Sq1.PRM.SMO*5));
    tc2=TCs;
    MSUA_resp=nan(SqTBU.PRM.nBoots,size(TCs,2));
    for n=1:SqTBU.PRM.nBoots
        for j=1:size(TCs,1)
            tc2(j,:)=circshift(TCs(j,:),round(rand(1,1)*size(TCs,2)));
        end
        MSUA_resp(n,:)=nanmean((tc2));
    end
    MSUA_act=nanmean((TCs));
    clearvars ZscoredMSUA
    for f=1:900
        mn_shf(f)=nanmean(MSUA_resp(:,f));
        sd_shf(f)=nanstd(MSUA_resp(:,f));
        temp=zscore([MSUA_resp(:,f) ; MSUA_act(f)]);
        ZscoredMSUA(f)=temp(end);
    end
    DumpMSUA(i,:)=MSUA_act;
    plot(MSUA_act,'color',CLR(i,:),'linewidth',2)
    hold on
    plot(mn_shf+BonFerri(900)*sd_shf,'color',[0.7 0.7 0.7])
    plot(mn_shf-BonFerri(900)*sd_shf,'color',[0.7 0.7 0.7])
    box off
    axis tight
    if i==1
        xlabel('Frame #')
        ylabel('MSUA (norm.)')
    end
    axis tight;%ylim([0.7 1.7])
    xticks([400 800]);
    xlim([1 900])
    MnRef=mn_shf;
    pkid=MSUA_act>mn_shf+BonFerri(900)*sd_shf;
    trid=MSUA_act<mn_shf-BonFerri(900)*sd_shf;
    SigPeaks(i,:)=nansum(pkid);
    SigTrs(i,:)=nansum(trid);
    PowPeaks(i,:)=nanmean(MSUA_act(pkid)- (mn_shf(pkid)+BonFerri(900)*sd_shf(pkid)) )/ nanmean(mn_shf)*100;
    PowTrs(i,:)=nanmean(-MSUA_act(trid)+  (mn_shf(trid)-BonFerri(900)*sd_shf(trid)) )/ nanmean(mn_shf)*100;
    
    XCorPeakHist(i,4,:)=circ_xcorr(MSUA_act,jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1),'corrcoef');
    [CorrCf(i,4,1),CorrCf(i,4,2)]=corrcoefwithNaNDiag(...
        MSUA_act,...
        jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1));
    
    %------------------------------Snake plot each spike gets equal vote-------------------------------------
    if i<4;subplot_CSP(4,4,15,6);else;subplot_CSP(4,4,16,7);end
    PYR=mintersect(IDAll{i},Broad2,Active2,TunedCells );
    TCs=SqTBU.FR_act(PYR,:)/delTCurrent;TCs=(colsmoothChinmay_Circ(TCs,Sq1.PRM.SMO*5));
    MSUA_act=nanmean((TCs));
    if i==1 || i==4
        yyaxis left
        plot(jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1),'color',[0.7 0.7 0.7])
    end
    yyaxis right
    hold on
    plot((MSUA_act),'color',CLR(i,:),'linewidth',1,'linestyle','-')
    box off
    XCorPeakHist(i,5,:)=circ_xcorr(MSUA_act,jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1),'corrcoef');
    [CorrCf(i,5,1),CorrCf(i,5,2)]=corrcoefwithNaNDiag(...
        MSUA_act,...
        jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1));
    ax=gca;ax.YAxis(1).Color = [0.7 0.7 0.7];ax.YAxis(2).Color = 'k';
    axis tight
    %---------------------------------------------------------------------
    figure(3)%Histogram of location of peaks, across the 7 brain regions
    set(gcf,'position',1e3*[ 0.3962    0.9786    1.5200    0.7248])
    
    %----------------- New stuff here about the F2F and the population overlap in adjecent frames
    if i==1;subplot_CSP(4,8,i,1);else;subplot(4,8,i);end
    PYR=mintersect(IDAll{i},Broad2,Active2,TunedCells );
    TCs=SqTBU.FR_act(PYR,:)/delTCurrent;
    TCs=(colsmoothChinmay_Circ(TCs,-Sq1.PRM.SMO));
    
    PopOverlap=nan(size(TCs,2),1);
    PopOverlap(1)=corrcoefwithNaNDiag(TCs(:,1),TCs(:,900));
    for ovr=2:900
        PopOverlap(ovr)=corrcoefwithNaNDiag(TCs(:,ovr-1),TCs(:,ovr));
    end
    PopOverlap=jmm_smooth_1d_cor(PopOverlap,Sq1.PRM.SMO,1);
    plot(jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1),'color',[0.7 0.7 0.7],'linewidth',2)
    if i==1;ylabel('F2F image correlation');end
    axis tight;ylim([0.65 1]);%yticks([0.75 1])
    yyaxis right
    hold on
    plot(PopOverlap,'color',CLR(i,:),'linewidth',1,'linestyle','-')
    box off
    if i==1;ylabel('Population overlap')
        xlabel('Frame #');end
    xticks([0 400 800]);
    ax=gca;ax.YAxis(1).Color = [0.7 0.7 0.7];ax.YAxis(2).Color = CLR(i,:);
    axis tight
    %ylim([0.994 1]);yticks([0.995 1])
    ylim([0.65 1]);yticks([0.75 1])
    
    XCorPeakHist(i,6,:)=circ_xcorr(PopOverlap,jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1),'corrcoef');
    [CorrCf(i,6,1),CorrCf(i,6,2)]=corrcoefwithNaNDiag(...
        PopOverlap,...
        jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1));
    DumpHistogram(i,1,:)=PopOverlap;
    title(TEXTsmall{i},'color',CLR(i,:));
    set(gca,'position',get(gca,'position')+[-(i-2)/1000 0 0.0 0 ]);
    %---------------------------------Peak location histogram--------------------
    if i==1;subplot_CSP(4,8,i+8,2);else;subplot(4,8,i+8);end
    yyaxis left
    plot(jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1),'color',[0.7 0.7 0.7],'linewidth',2)
    if i==1;ylabel('F2F image correlation');end
    axis tight;ylim([0.65 1]);%yticks([0.75 1])
    yyaxis right
    [a,b]=histperc(DumpPeakProp{i,12},BinFrame);
    XCorPeakHist(i,1,:)=circ_xcorr(interp1(center(BinFrame),jmm_smooth_1d_cor(a,2,1),center(Sq1.PRM.FrameBins)),jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1),'corrcoef');
    [CorrCf(i,1,1),CorrCf(i,1,2)]=corrcoefwithNaNDiag(...
        interp1(center(BinFrame),jmm_smooth_1d_cor(a,2,1),center(Sq1.PRM.FrameBins)),...
        jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1));
    [~ ,p_chisq(i,1)]= chi2gof(DumpPeakProp{i,12},'CDF',pd);
    stairs(b,jmm_smooth_1d_cor(a,2,1),'linewidth',2,'color',CLR(i,:));hold on
    box off; axis tight; ylim([0 5]);
    xticks([0 400 800])
    if i==1;xlabel('Frame #');ylabel('% of peaks');
    end
    ax=gca;ax.YAxis(1).Color = [0.7 0.7 0.7];ax.YAxis(2).Color = CLR(i,:);
    set(gca,'position',get(gca,'position')+[-(i-2)/1000 0 0.0 0 ]);
    DumpHistogram(i,2,:)=interp1(center(BinFrame),jmm_smooth_1d_cor(a,2,1),center(Sq1.PRM.FrameBins));
    %-----------------------------Peak width histogram--------------------
    if i==1;subplot_CSP(4,8,i+16,3);else;subplot(4,8,i+16);end
    yyaxis left
    plot(jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1),'color',[0.7 0.7 0.7],'linewidth',2)
    if i==1;ylabel('F2F image correlation');end
    axis tight;ylim([0.65 1]);%yticks([0.75 1]);
    yyaxis right
    [a,b]=depHist(DumpPeakProp{i,5},DumpPeakProp{i,12},BinFrame,'Given',0,@nanmedian,0,'none',1);
    [a1]=depHist(DumpPeakProp{i,5},DumpPeakProp{i,12},Sq1.PRM.FrameBins(1:5:end),'Given',0,@nanmax,0,'none',1);
    [a2]=depHist(DumpPeakProp{i,5},DumpPeakProp{i,12},Sq1.PRM.FrameBins(1:5:end),'Given',0,@nanmin,0,'none',1);
    DumpMegaScaleRatio(i,:)=a1./a2;
    XCorPeakHist(i,2,:)=circ_xcorr(interp1(center(BinFrame),jmm_smooth_1d_cor(a,2,1),center(Sq1.PRM.FrameBins)),jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1),'corrcoef');
    [CorrCf(i,2,1),CorrCf(i,2,2)]=corrcoefwithNaNDiag(...
        interp1(center(BinFrame),jmm_smooth_1d_cor(a,2,1),center(Sq1.PRM.FrameBins)),...
        jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1));
    stairs(b,jmm_smooth_1d_cor(a,2,1),'linewidth',2,'color',CLR(i,:));hold on
    box off; axis tight;
    
    if i==1;xlabel('Frame #');ylabel('Median width (sec)');end
    ax=gca;ax.YAxis(1).Color = [0.7 0.7 0.7];ax.YAxis(2).Color = CLR(i,:);
    yyaxis right;ylim([0.1 12]);set(gca,'yscale','log');
    yticks([0.1 1 10])
    [~ ,p_chisq(i,2)]= chi2gof(DumpPeakProp{i,5},'CDF',pd);
    xticks([0 400 800])
    set(gca,'position',get(gca,'position')+[-(i-2)/1000 0 0.0 0 ]);
    DumpHistogram(i,3,:)=interp1(center(BinFrame),jmm_smooth_1d_cor(a,2,1),center(Sq1.PRM.FrameBins));
    %-----------------------------------All units--------------------------------
    if i==1;subplot_CSP(4,8,i+24,4);else;subplot(4,8,i+24);end
    AllSpikes=mintersect(IDAll{i});
    TCs=SqTBU.FR_act(AllSpikes,:)/delTCurrent;TCs=(colsmoothChinmay_Circ(TCs,Sq1.PRM.SMO*5));
    tc2=TCs;
    MSUA_resp=nan(SqTBU.PRM.nBoots,size(TCs,2));
    for n=1:SqTBU.PRM.nBoots
        for j=1:size(TCs,1)
            tc2(j,:)=circshift(TCs(j,:),round(rand(1,1)*size(TCs,2)));
        end
        MSUA_resp(n,:)=nanmean((tc2));
    end
    MSUA_act=nanmean((TCs));
    clearvars ZscoredMSUA
    for f=1:900
        mn_shf(f)=nanmean(MSUA_resp(:,f));
        sd_shf(f)=nanstd(MSUA_resp(:,f));
        temp=zscore([MSUA_resp(:,f) ; MSUA_act(f)]);
        ZscoredMSUA(f)=temp(end);
    end
    %DumpMSUA(i,:)=ZscoredMSUA;
    plot(MSUA_act,'color',CLR(i,:),'linewidth',2)
    hold on
    plot(mn_shf+BonFerri(900)*sd_shf,'color',[0.7 0.7 0.7],'linewidth',2)
    plot(mn_shf-BonFerri(900)*sd_shf,'color',[0.7 0.7 0.7],'linewidth',2)
    box off
    axis tight
    XCorPeakHist(i,3,:)=circ_xcorr(MSUA_act,jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1),'corrcoef');
    [CorrCf(i,3,1),CorrCf(i,3,2)]=corrcoefwithNaNDiag(...
        MSUA_act,...
        jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1));
    xticks([0 400 800])
    if i==1;ylabel('Firing rate (Hz)')
        xlabel('Frame #');end
    DumpHistogram(i,4,:)=MSUA_act;
    
    MnRef=mn_shf;
    pkid=MSUA_act>mn_shf+BonFerri(900)*sd_shf;
    trid=MSUA_act<mn_shf-BonFerri(900)*sd_shf;
    SigPeaks_allFR(i,:)=nansum(pkid);
    SigTrs_allFR(i,:)=nansum(trid);
    PowPeaks_allFR(i,:)=nanmean(MSUA_act(pkid)- (mn_shf(pkid)+BonFerri(900)*sd_shf(pkid)) )/ nanmean(mn_shf)*100;
    PowTrs_allFR(i,:)=nanmean(-MSUA_act(trid)+  (mn_shf(trid)-BonFerri(900)*sd_shf(trid)) )/ nanmean(mn_shf)*100;
    set(gca,'position',get(gca,'position')+[-(i-2)/1000 0 0.0 0 ]);
    
end
figure(1)
%p=subplot(8,2,15);
subplot_CSP(4,4,13,4);
b=bar(1:7,[SigPeaks SigTrs]);
b(1).FaceColor=CLRPeak;b(1).EdgeColor=CLRPeak;
b(2).FaceColor=CLRTrough;b(2).EdgeColor=CLRTrough;
xticks([1:7])
xticklabels(TEXTsmall)
ylabel('# Frames with significant deviation')
box off
L=legend('Above chance','Below');L.ItemTokenSize(1)=7;legend boxoff;box off;
%title('MSUA extrema are rare in hippocampal regions, so dont read much into that')
xtickangle(25)

subplot_CSP(4,4,14,5);
%p=subplot(5,2,10);
b=bar(1:7,[PowPeaks PowTrs]);
b(1).FaceColor=CLRPeak;b(1).EdgeColor=CLRPeak;
b(2).FaceColor=CLRTrough;b(2).EdgeColor=CLRTrough;
xticks([1:7])
xticklabels(TEXTsmall)
ylabel('Firing deviation (%)')
box off
xtickangle(25)
%title('MSUA extrema are rare in hippocampal regions, so dont read much into that')

figure(3)
clearvars Corrmat PValMat
for k=1:4
for i=1:7
for j=(i):7

   [tempcrr(1,:,:),tempcrr(2,:,:)]= PartialCorrwithNan([squeeze(DumpHistogram(i,k,:)),squeeze(DumpHistogram(j,k,:)),...
    jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1)']);
   [tempcr(1),tempcr(2)]= corrcoefwithNaNDiag(squeeze(DumpHistogram(i,k,:)),squeeze(DumpHistogram(j,k,:)));
   [tempcr2(1),tempcr2(2)]= corrcoefwithNaNDiag(squeeze(DumpHistogram(i,k,:)),jmm_smooth_1d_cor(MV.Coh,Sq1.PRM.SMO,1)');
Corrmat(k,i,j)=tempcrr(1,1,2);PValMat(k,i,j)=tempcrr(2,1,2);
Corrmat(k,j,i)=tempcr(1);PValMat(k,j,i)=tempcr(2);
Corrmat(k,i,i)=tempcr2(1);PValMat(k,i,i)=tempcr2(2);
end
end
subplot(4,8,8*k);
imagescwithnan(1:7,1:7,squeeze(Corrmat(k,:,:)),jet,[ 0 0 0]);if k==1;colorbar('westOutside');end;caxis([-1 1])
xticks([1:7]);yticks([1:7]);
xticklabels([TEXTsmall]);yticklabels([TEXTsmall]);xtickangle(45)
if k==1;title('Correlation matrix');end
end
%% MSUA ED for scrambled data
figure('Position',[507.4000   38.6000  992.8000  695.2000])
clc
SqTBU=Sc1;
clearvars PowPeaks PowTrs SigPeaks SigTrs r s AvgNonDiag
for i=1:3
    p=subplot(5,3,[i i+3]);
    p1=get(p,'position');
    PYR=mintersect(find(SqTBU.Zspa>2),IDAll{i},Broad2,Active2 );
    TCs=SqTBU.FR_act(PYR,:);
    TCs=colmaxNormalize(colsmoothChinmay_Circ(TCs,Sq1.PRM.SMO*5));
    temp=TCs(:);
    AvgNonDiag(i)=nanmean(temp(temp~=1));
    SnakePlotMaker_removeNan(center(SqTBU.PRM.FrameBins),TCs);
    if i==1;xlabel('Frame #');ylabel('Cell #');
        annotation('textbox',[p1(1)-p1(3)*0.3  p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],...
            'string','a','EdgeColor','none','FontSize',18,'FontWeight','bold')
    end
    colorbar off;axis tight;ylim([0 inf])
    xticks([400 800]);%yticks([0 floor(size(TCs(~isnan(TCs(:,1)),:),1)/10)*10])
    box off
    title(TEXTsmall{i},'color',CLR(i,:))
    ax=gca;ax.YAxis(1).Color = CLR(i,:);ax.XAxis(1).Color = CLR(i,:);
    
    p=subplot(5,3,[i+6 i+9]);
    p1=get(p,'position');
    TCs=SqTBU.FR_act(PYR,:);
    %TCs=colmeanNormalize(colsmoothChinmay_Circ(TCs,Sq1.PRM.SMO*5));
    TCs=colmaxNormalize(colsmoothChinmay_Circ(TCs,Sq1.PRM.SMO*5));
    tc2=TCs;
    MSUA_resp=nan(SqTBU.PRM.nBoots,size(TCs,2));
    for n=1:SqTBU.PRM.nBoots
        for j=1:size(TCs,1)
            tc2(j,:)=circshift(TCs(j,:),round(rand(1,1)*size(TCs,2)));
        end
        MSUA_resp(n,:)=nanmean((tc2));
    end
    MSUA_act=nanmean((TCs));
    clearvars ZscoredMSUA
    for f=1:900
        mn_shf(f)=nanmean(MSUA_resp(:,f));
        sd_shf(f)=nanstd(MSUA_resp(:,f));
        temp=zscore([MSUA_resp(:,f) ; MSUA_act(f)]);
        ZscoredMSUA(f)=temp(end);
    end
    DumpMSUA(i,:)=MSUA_act;
    plot(MSUA_act,'color',CLR(i,:),'linewidth',2)
    hold on
    plot(mn_shf+BonFerri(900)*sd_shf,'color',[0.7 0.7 0.7])
    plot(mn_shf-BonFerri(900)*sd_shf,'color',[0.7 0.7 0.7])
    box off
    axis tight
    if i==1
        xlabel('Frame #')
        ylabel('MSUA (norm.)')
    end
    %ylim([0.45 0.7])
    xticks([400 800]);
    xlim([1 900])
    MnRef=mn_shf;
    pkid=MSUA_act>mn_shf+BonFerri(900)*sd_shf;
    trid=MSUA_act<mn_shf-BonFerri(900)*sd_shf;
    SigPeaks(i,:)=nansum(pkid);
    SigTrs(i,:)=nansum(trid);
    PowPeaks(i,:)=nanmean(MSUA_act(pkid)- (mn_shf(pkid)+BonFerri(900)*sd_shf(pkid)) )/ nanmean(mn_shf)*100;
    PowTrs(i,:)=nanmean(-MSUA_act(trid)+  (mn_shf(trid)-BonFerri(900)*sd_shf(trid)) )/ nanmean(mn_shf)*100;
    
    XCorPeakHist_sca(i,4,:)=circ_xcorr(MSUA_act,jmm_smooth_1d_cor(Coh2,Sq1.PRM.SMO,1),'corrcoef');
    [r_sca(i),s_sca(i)]=corrcoefwithNaNDiag(...
        MSUA_act,...
        jmm_smooth_1d_cor(Coh2,-Sq1.PRM.SMO,1));
end

p=subplot(5,3,13);
p1=get(p,'position');
annotation('textbox',[p1(1)-p1(3)*0.3  p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],...
    'string','c','EdgeColor','none','FontSize',18,'FontWeight','bold')
b=bar(1:3,[SigPeaks SigTrs]);
b(1).FaceColor=CLRPeak;b(1).EdgeColor=CLRPeak;
b(2).FaceColor=CLRTrough;b(2).EdgeColor=CLRTrough;
xticks([1:3])
xticklabels(TEXTsmall(1:3))
ylabel('# Frames with significant deviation')
box off
L=legend('Above chance','Below');L.ItemTokenSize(1)=7;legend boxoff;box off;
%title('MSUA extrema are rare in hippocampal regions, so dont read much into that')
xtickangle(25)

p=subplot(5,3,14);
p1=get(p,'position');
annotation('textbox',[p1(1)-p1(3)*0.3  p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],...
    'string','d','EdgeColor','none','FontSize',18,'FontWeight','bold')
b=bar(1:3,[PowPeaks PowTrs]);
b(1).FaceColor=CLRPeak;b(1).EdgeColor=CLRPeak;
b(2).FaceColor=CLRTrough;b(2).EdgeColor=CLRTrough;
xticks([1:3])
xticklabels(TEXTsmall(1:3))
ylabel('Firing deviation (%)')
box off
xtickangle(25)

p=subplot(5,3,15);
p1=get(p,'position');
annotation('textbox',[p1(1)-p1(3)*0.3  p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],...
    'string','e','EdgeColor','none','FontSize',18,'FontWeight','bold')
plot(1:900,jmm_smooth_1d_cor(Coh2,-Sq1.PRM.SMO,1),'color',[[0.5 0.5 1]])
box off
xlabel('Frame #')
ylabel('Frame to frame correlation')
axis tight

%% F4 Check the amount of tuning for sequential and scrambled movies
figure('position',1e3*[0.3850    0.1074    1.0128    0.7000])
clearvars storeAll storeMiddle storeSca ks_pvalue ks_pvalue_wdShf storeScaCompare
load('Movie_Scra_Corrcoef_May1_2022.mat');
for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq1.Zspa)),Active2 );
    storeAll(i,1:3)=[nanlength(PYR) nansum(Sq1.Zspa(PYR)>2) FracAbove2(Sq1.Zspa(PYR))] ;
    PYR2d=mintersect(Broad2,IDAll{i},find(~isnan(Sq3.Zspa)),Active2 );
    storeMiddle(i,1:3)=[nanlength(PYR2d) nansum(Sq3.Zspa(PYR2d)>2) FracAbove2(Sq3.Zspa(PYR2d))] ;
    PYR2d=mintersect(Broad2,IDAll{i},find(~isnan(Sc1.Zspa)),Active2 );
    storeSca(i,1:3)=[nanlength(PYR2d) nansum(Sc1.Zspa(PYR2d)>2) FracAbove2(Sc1.Zspa(PYR2d))] ;
    storeScaCompare(i,1:3)=[nanlength(PYR2d) nansum(Sc1.Zspa_shf(PYR2d)>2) FracAbove2(Sc1.Zspa_shf(PYR2d))] ;
    [~,ks_pvalue(i)]=kstest2(Sq3.Zspa(PYR),Sc1.Zspa(PYR2d));
    [~,ks_pvalue_wdShf(i)]=kstest2(Sc1.Zspa_shf(PYR2d),Sc1.Zspa(PYR2d));
end

subplot_CSP(4,3,4);
b=bar(1:3,[storeAll(1:3,end) storeMiddle(1:3,end) storeSca(1:3,end)]*100);
box off;
b(1).FaceColor=ColorSeq;b(1).EdgeColor=b(1).FaceColor;
b(2).FaceColor=ColorSeq2;b(2).EdgeColor=b(2).FaceColor;
b(3).FaceColor=ColorSca;b(3).EdgeColor=b(3).FaceColor;
ylabel('Visual tuning(%)')
line([1 3],[2.25 2.25],'color',[0 0 0],'linestyle','--')
ylim([0 100])
yyaxis right
b=bar(4:7,[storeAll(4:7,end) storeMiddle(4:7,end) storeSca(4:7,end)]*100);
line([4 7],[2.25 2.25],'color',[0 0 0],'linestyle','--')
b(1).FaceColor=ColorSeq;b(1).EdgeColor=b(1).FaceColor;
b(2).FaceColor=ColorSeq2;b(2).EdgeColor=b(2).FaceColor;
b(3).FaceColor=ColorSca;b(3).EdgeColor=b(3).FaceColor;
xticks([1:7])
xticklabels([TEXTsmall])
xtickangle(20)
ax=gca;ax.YAxis(1).Color ='k';ax.YAxis(2).Color ='k';
line([3.5 3.5],ylim,'linestyle','-.')
box off
ylabel('Hippocampal tuning (%)')
L=legend('Continuous (all)','Continuous (subset)','Scrambled','Chance level');
legend boxoff;L.ItemTokenSize(1)=8;
ylim([0 70])

subplot_CSP(4,3,1);
if ~exist( 'Coh1','var')  && ~exist( 'Coh2','var')
    for i=1:900
        if i>1
            Coh1(i)=corrcoefwithNaNDiag(MV.Movie(:,:,i-1),MV.Movie(:,:,i));
            Coh2(i)=corrcoefwithNaNDiag(MV.Movie_Shf(:,:,i-1),MV.Movie_Shf(:,:,i));
        else
            Coh1(i)=corrcoefwithNaNDiag(MV.Movie(:,:,900),MV.Movie(:,:,i));
            Coh2(i)=corrcoefwithNaNDiag(MV.Movie_Shf(:,:,900),MV.Movie_Shf(:,:,i));
        end
    end
end
plot(jmm_smooth_1d_cor(Coh1,-Sq1.PRM.SMO,1),'color',ColorSeq2,'linewidth',1);hold on;
plot(jmm_smooth_1d_cor(Coh2,-Sq1.PRM.SMO,1),'color',ColorSca,'linewidth',1)
axis tight;box off;ylabel('F2F correlation');xlabel('Frame #')
L=legend('Continuous','Scrambled');legend boxoff;L.ItemTokenSize(1)=8;
ylim([-0.4 1.05])

%MSUA response and population overlap also go in this figure. Also show
%example cells 1 each brain region
A2sca=[A2(4,2); A2(6,1)];clrID=[4 6];
SPLoc=[2 5;3 6];
UseAllData=1;
for j=1:size(A2sca,1)
    for i=1:size(A2sca,2)
        idcurrent=A2sca(j,i);
        ses=DumpNS.ExpID(idcurrent);
        FileName=[D(ses).folder '\' D(ses).name];
        UnitID=DumpNS.UnitID(idcurrent);
        FR_actC=Sc1.FR_act(idcurrent,:)/delTCurrent;
        FR_shf_meanC=Sc1.FR_shf_mean(idcurrent,:)/delTCurrent;
        FR_shf_stdC=Sc1.FR_shf_std(idcurrent,:)/delTCurrent;
        if i==1;s1=subplot(6,3,SPLoc(j,:));else;s1=subplot(6,3,SPLoc(j,:)+6);end
        LocIN=[s1.Position];if i==2; LocIN= LocIN+[0 0.0 0 0];end
        if i==1 && j==1
            LabelBool=1;
            annotation('textbox',[LocIN(1)-LocIN(3)*0.20 LocIN(2)+LocIN(4)*1.11 LocIN(3)/10 LocIN(4)/10],...
                'string','d','EdgeColor','none','FontSize',18,'FontWeight','bold')
        else
            LabelBool=0;
        end
        l1= Movie_PlotScatter_PSTH_psychedelic_shuffled(UnitID,FileName,Sc1.PRM.TCutoff,Sc1.PRM.RVCut,Sc1.PRM.FrameBins,...
            FR_actC,FR_shf_meanC,FR_shf_stdC,LocIN,LabelBool,CLR(clrID(j),:),UseAllData);
        set(gcf,'CurrentAxes',l1)
        if i==1;title(['' DumpNS.Region{idcurrent}],'color',CLR(clrID(j),:));end
        %title([num2str(idcurrent) '-' DumpNS.Region(idcurrent)]);
    end
    % annotation('rectangle',[s1(1)-s1(3)*0.2 s1(2)+s1(4)*0.2 s1(3)/10 s1(4)/10],'Color',CLR(j+3,:),'Linewidth',2)
end
annotation('rectangle',[0.33 0.66 0.66 0.33],'Color',ColorSca,'linewidth',2)

subplot_CSP(4,3,[7 10],3);
clearvars SummaryTuningChange DumpTuningChange
for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq3.Zspa) & ~isnan(Sc1.Zspa) & (Sq3.Zspa>2 | Sc1.Zspa>2)),Active2,Active2_sca );
    MetricCurrent=100*(Sq3.Zspa(PYR) - Sc1.Zspa(PYR))./max(Sq3.Zspa(PYR),Sc1.Zspa(PYR));
    [a,b]=ecdf( MetricCurrent );plot(b,a,'color',CLR(i,:),'linewidth',2)
    
    DumpTuningChange{i}=MetricCurrent ;
    SummaryTuningChange(i,:)=[nanmedian(MetricCurrent) nanmean(MetricCurrent) ...
        nanstd(MetricCurrent) nanstd(MetricCurrent)/sqrt(nanlength(MetricCurrent))];
    hold on
end
xlim([-150 200]);
L=legend(TEXTsmall) ;legend boxoff;L.ItemTokenSize(1)=8;
xlabel('Tuning change = 100*(Z_{continuous} - Z_{scrambled})/max(Z_{continuous},Z_{scrambled}) (%)') ;%100*(Z_{sequential} - Z_{scrambled})/max(Z_{sequential},Z_{scrambled})
ylabel('Probability')
box off
grid on

%{
 %__________________________________ Decoding the movie and scramble
TRW_seq=load('Movie_FRStore_TrWise_AllTr_Seq.mat');
TRW_sca=load('Movie_FRStore_TrWise_AllTr_Sca.mat');
TRW_sca_rearr=TRW_sca;
id=cell2mat(IDAll([1:2]));TRW_sca_rearr.TCDump_TrWise(id,:,[899 900 1:898 ])=TRW_sca.TCDump_TrWise(id,:,1:900);TRW_sca_rearr.TCDump_TrWise(id,:,MV.ShuffledFrameIDs)=TRW_sca_rearr.TCDump_TrWise(id,:,:);
id=IDAll{3};TRW_sca_rearr.TCDump_TrWise(id,:,[898 899 900 1:897 ])=TRW_sca_rearr.TCDump_TrWise(id,:,1:900);TRW_sca_rearr.TCDump_TrWise(id,:,MV.ShuffledFrameIDs)=TRW_sca_rearr.TCDump_TrWise(id,:,:);

Error_Act=nan(3,7,60,3,2);%Error_Bayes_Act=nan(7,1);
Error_Shf=nan(3,7,60,3,2);%Error_Bayes_Shf=nan(7,nBoots);
Error_Act2=nan(3,7,60,3,2);%Error_Bayes_Act=nan(7,1);
Error_Shf2=nan(3,7,60,3,2);%Error_Bayes_Shf=nan(7,nBoots);
Out_Act=nan(3,7,60,SizeBins,3,2);%Out_Bayes_Act=nan(7,SizeBins);
Out_Shf=nan(3,7,60,SizeBins,3,2);%Out_Bayes_shf=nan(7,nBoots,SizeBins);
ProbOut=nan(7,SizeBins,SizeBins);
nCellsDecoding=nan(3,7,60);
MinCellCount=150;
downsampleBool=[1 0];
SMOLevels=[2 30];
for i=1:7
    for Type=2:3
        if Type==1
            TRW_touse=TRW_seq;
            TRAll=1:60;
            idc=mintersect(Broad2,IDAll{i},find(Sq1.Zspa>2),Active2);
        end
        if Type==2
            TRW_touse=TRW_seq;
            TRAll=21:40;
            idc=mintersect(Broad2,IDAll{i},find(Sq3.Zspa>2),Active2);
        end
        if Type==3
            TRW_touse=TRW_sca;
            TRAll=1:20;
            idc=mintersect(Broad2,IDAll{i},find(Sc1.Zspa>2),Active2_sca);
        end
        for DownSample=[2 ]%1:2%2
            for smo=1:2
                SMOLevel=SMOLevels(smo);
                
                for tr=TRAll
                    [i tr]
                    clearvars b b2
                    ENC=squeeze(nanmean( squeeze(TRW_touse.TCDump_TrWise(idc,setdiff(TRAll,tr),:)) , 2));ENC=colNormalize01(colsmoothChinmay_Circ(ENC,SMOLevel));
                    DEC=squeeze(( squeeze(TRW_touse.TCDump_TrWise(idc,tr,:))));DEC=colNormalize01(colsmoothChinmay_Circ(DEC,SMOLevel));
                    temp=~isnan(ENC(:,1)) & ~isnan(DEC(:,1));
                    nCellsDecoding(Type,i,tr)=sum(temp);
                    if downsampleBool(DownSample)==1
                        allcells=find(temp==1);chosenCells=allcells(randperm(length(allcells),MinCellCount));
                    else
                        chosenCells=find(temp==1);
                    end
                    ENC2=ENC(chosenCells,:);DEC2=DEC(chosenCells,:);
                    Out2=VecCorrCoef_v3_fast(ENC2,DEC2);Out=Out2;Out(Out2<0)=0;
                    [~,b]=nanmax(Out');
                    Error_Act(Type,i,tr,smo,DownSample)=sqrt(sum( (b-(1:900)).^2 ))/SizeBins;
                    Error_Act2(Type,i,tr,smo,DownSample)=nanmean(abs(b-(1:900)));
                    Out_Act(Type,i,tr,:,smo,DownSample)=b;
                    for n=1:1
                        clearvars b b2
                        RNDPRM=randperm(size(ENC2,1));
                        ENC2s=ENC2(RNDPRM,:);
                        Out2=VecCorrCoef_v3_fast(ENC2s,DEC2);Out=Out2;Out(Out2<0)=0;
                        [~,b]=nanmax(Out');
                        Error_Shf(Type,i,tr,smo,DownSample)=sqrt(sum( (b-(1:900)).^2 ))/SizeBins;
                        Error_Shf2(Type,i,tr,smo,DownSample)=nanmean(abs(b-(1:900)));
                        Out_Shf(Type,i,tr,:,smo,DownSample)=b;
                    end
                end
            end
        end
    end
end
save('Movie_decodingData_Feb9_2023_v3.mat','Error_Act','Error_Act2','Error_Shf','Error_Shf2','Out_Act','Out_Shf','SMOLevel','downsampleBool','nCellsDecoding','-v7.3');
%}
load('Movie_decodingData_Feb9_2023_v3.mat');
subplot_CSP(5,3,8);
h1=BoxPlot_Chinmay(1:7,(squeeze(Error_Act2(1,:,:,1,2)))',ColorSeq,1);
hold on
h2=BoxPlot_Chinmay(1:7,(squeeze(Error_Shf2(1,:,:,1,2)))',[0.7 0.7 0.7],1);
hold on
h3=plot(1:7,nanmean((squeeze(Error_Act2(1,:,:,1,1)))'),'*','color',[1 1 0],'linewidth',3);
box off;axis tight;ylim([0 400])
xticks(1:7);xticklabels(TEXTsmall);xtickangle(20)
ylabel('Decoding error (# frames)')
legend([h1(1) h2(1) h3],{'Continuous (all trials, all cells)','Shuffle','Continuous (all trials, 150 cells'})
L=legend;legend boxoff;L.ItemTokenSize(1)=8;

subplot_CSP(5,3,9);
h1=BoxPlot_Chinmay(1:7,(squeeze(Error_Act2(2,:,:,1,2)))',ColorSeq2,1);
hold on
h2=BoxPlot_Chinmay(1:7,(squeeze(Error_Act2(3,:,:,1,2)))',ColorSca,1);
hold on
h3=plot(1:7,nanmean((squeeze(Error_Act2(1,:,:,1,1)))'),'*','color',[1 1 0],'linewidth',3);
box off;axis tight;ylim([0 400])
xticks(1:7);xticklabels(TEXTsmall);xtickangle(20)
legend([h1(1) h2(1)],{'Continuous (20 trials)','Scrambled'})
L=legend;legend boxoff;L.ItemTokenSize(1)=8;

for i=1:7
    [~,rnk_p1(i)]=kstest2(squeeze(Error_Act2(1,i,:,1,2))',squeeze(Error_Shf2(1,i,:,1,2))');
    [~,rnk_p2(i)]=kstest2(squeeze(Error_Act2(2,i,:,1,2))',squeeze(Error_Act2(3,i,:,1,2))');
    [~,rnk_p3(i)]=kstest2(squeeze(Error_Shf2(3,i,:,1,2))',squeeze(Error_Act2(3,i,:,1,2))');
    [~,rnk_p4(i)]=kstest2(squeeze(Error_Act2(1,i,:,1,1))',squeeze(Error_Shf2(1,i,:,1,1))');
end

%-------------- Now I want to plot the EG cells which could and could not be rearranged
SmoLevel=Sq1.PRM.SMO*Sq5_hr.PRM.Upsample;
XBins=center(Sq5_hr_middleTrials.PRM.FrameBins);
ChosenID=[14314 50222];%ca1 13913, 14177, 17217, 49384
CLRID=[2 6];
LocID=[11 12; 14 15];
for i=1:2
    tc1=Sq5_hr_middleTrials.SpCnt_actual((ChosenID(i)),:)/SpCnt2FRScaler;
    tc2=Sc5_hr.SpCnt_actual((ChosenID(i)),:)/SpCnt2FRScaler;
    [bestCRCF,bestLatency]=max(CRCF_dump_hr(ChosenID(i),:));
    tc1_shifted=circshift(tc1,Timeshift(bestLatency));
    tc2_shifted=circshift(tc2,Timeshift(bestLatency));
    tc2_reassigned(ShuffledFrameIDs_hr)=tc2_shifted;
    p=subplot(5,3,LocID(i,:));
    plot(XBins,jmm_smooth_1d_cor_circ(tc2_shifted,SmoLevel,1),'color',[0.8 0.8 0.8],'linewidth',1);
    hold on;
    plot(XBins,jmm_smooth_1d_cor_circ(tc1_shifted,SmoLevel,1),'color',ColorSeq2,'linewidth',1);
    plot(XBins,jmm_smooth_1d_cor_circ(tc2_reassigned,SmoLevel,1),'color',ColorSca,'linewidth',1);
    box off
    axis tight
    xticks([400 800])
    if i==1;ylabel('FR (sp/sec)');
        L=legend('Scrambled - as is','Continuous','Scrambled - rearranged');
        legend boxoff;L.ItemTokenSize(1)=8;
        p1=get(p,'position');
        annotation('textbox',[p1(1)-p1(3)/10 p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],'string','g',...
            'EdgeColor','none','FontSize',18,'FontWeight','bold')
        xlabel('Frame #')
    end
    title([TEXTsmall{CLRID(i)} ' r=' num2str(roundsd(bestCRCF,2)) ],'color',CLR(CLRID(i),:));
    ax=gca;ax.YAxis.Color = CLR(CLRID(i),:);ax.XAxis.Color = CLR(CLRID(i),:);
end

%% Create the effective tuning curves in Even and odd trials, starting from the TrialWise FR matrices and create that plot

TRW_seq=load('Movie_FRStore_TrWise_AllTr_Seq.mat');
TRW_sca=load('Movie_FRStore_TrWise_AllTr_Sca.mat');

figure('Position',1e3*[0.1026    0.1210    1.2744    0.6320])
VecDiags=-899:899;
for i=1:7
    for Switch=1:2
        idc=mintersect(Broad2,IDAll{i},Active2,TunedCells);
        if Switch==1
            idc=mintersect(Broad2,IDAll{i},Active2,TunedCells);
            %TCShf=TCDump_TrWise(idc,:,:);
        else
            %TCShf=Allen_Movie_TrWiseFRMat_CircShifter(TCDump_TrWise(idc,:,:));
            idc=mintersect(Broad2,IDAll{i},Active2,find(Sq1.Zspa<2));
        end
        
        TC1=colmeanNormalize(colsmoothChinmay_Circ(squeeze(nanmean(TRW_seq.TCDump_TrWise(idc,2:2:end,:),2)),Sq1.PRM.SMO));
        TC2=colmeanNormalize(colsmoothChinmay_Circ(squeeze(nanmean(TRW_seq.TCDump_TrWise(idc,1:2:end,:),2)),Sq1.PRM.SMO));
        TC1(isnan(TC1))=0;TC2(isnan(TC2))=0;
        if Switch==1
            p=subplot(6,7,[i i+7]);
            p1=get(p,'position');
            if i==1;
                annotation('textbox',[p1(1)-p1(3)*0.3 p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],...
                    'string','a','EdgeColor','none','FontSize',18,'FontWeight','bold')
            end
            O=VecCorrCoef_v3_fast(TC1,TC2);
            imagescwithnan(center(Sq1.PRM.FrameBins),center(Sq1.PRM.FrameBins),(O'),jet,[1 1 1]);
            ax=gca;ax.XAxis.Color = CLR(i,:);ax.YAxis.Color = CLR(i,:);
            tempo=caxis;
            title([TEXTsmall{i} ' - Tuned'],'color',CLR(i,:))
        else
            p=subplot(6,7,[i+14 i+21]);
            p1=get(p,'position');
            if i==1;
                annotation('textbox',[p1(1)-p1(3)*0.3 p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],...
                    'string','b','EdgeColor','none','FontSize',18,'FontWeight','bold')
            end
            O=VecCorrCoef_v3_fast(TC1,TC2);
            imagescwithnan(center(Sq1.PRM.FrameBins),center(Sq1.PRM.FrameBins),(O'),jet,[1 1 1]);
            ax=gca;ax.XAxis.Color = [0.7 0.7 0.7];ax.YAxis.Color = [0.7 0.7 0.7];
            %caxis([tempo])
            title([TEXTsmall{i} ' - Untuned'],'color',[0.7 0.7 0.7])
            
        end
        %colorbar('SouthOutside')
        %caxis([-1 1])
        %if i<4;    caxis([0.36 1]);    else;        caxis([0.6 1.05]);    end
        axis square
        if i==1
            line([1 600],[301 900],'color','k','linestyle','--','linewidth',2);
            line([301 900],[1 600],'color','k','linestyle','--','linewidth',2);
            if Switch==1
                cbr=colorbar('SouthOutside');
                ylabel(cbr,'Pop. overlap','FontSize',12);
            end
        end
        caxis([-0.4 1])
        %Type 2 The various diagonals
        clearvars Diags1 Diags2 Diags3
        for k=-899:899
            Diags1(k+900)=nanmean(diag(O,k));
            Diags2(k+900)=nanstd(diag(O,k));
            Diags3(k+900)=nanlength(diag(O,k));
        end
        %figure(2)
        p=subplot(6,7,[i+28 i+35]);
        p1=get(p,'position');
        if i==1;
            annotation('textbox',[p1(1)-p1(3)*0.3 p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],...
                'string','c','EdgeColor','none','FontSize',18,'FontWeight','bold')
        end
        if Switch==1
            shadedErrorBar(VecDiags,Diags1,Diags2./sqrt(Diags3),CLR(i,:));
            hold on
        else
            shadedErrorBar(VecDiags,Diags1,Diags2./sqrt(Diags3),[0.7 0.7 0.7]);
            hold on
        end
        axis tight;box off;ylim([-0.2 0.9])
        if i==1;
            xlabel('Diagonal number');ylabel('Avg. correlation');
            line([-300 -300],ylim,'color','k','linestyle','--','linewidth',2);
            line([300 300],ylim,'color','k','linestyle','--','linewidth',2);
        end
        if Switch==1
            DiagsOut(i,1,:)=Diags1;
        else
            DiagsOut_shf(i,2,:)=Diags1;
        end
        xlim([-450 450])
        %if i<4;ylim([0.5 1.25]);else;ylim([0.85 1.15]);end
        
    end
end

%%% Scrambled movie population vector overlap

figure('Position',1e3*[0.1026    0.1210    1.2744    0.6320])
VecDiags=-899:899;
for i=1:7
    for Switch=1:2
        
        if Switch==1
            idc=mintersect(Broad2,IDAll{i},Active2_sca,find(Sc1.Zspa>2));
            %TCShf=TCDump_TrWise(idc,:,:);
        else
            %TCShf=Allen_Movie_TrWiseFRMat_CircShifter(TCDump_TrWise(idc,:,:));
            idc=mintersect(Broad2,IDAll{i},Active2_sca,find(Sc1.Zspa<2));
        end
        
        TC1=colmeanNormalize(colsmoothChinmay_Circ(squeeze(nanmean(TRW_sca.TCDump_TrWise(idc,2:2:end,:),2)),Sq1.PRM.SMO));
        TC2=colmeanNormalize(colsmoothChinmay_Circ(squeeze(nanmean(TRW_sca.TCDump_TrWise(idc,1:2:end,:),2)),Sq1.PRM.SMO));
        TC1(isnan(TC1))=0;TC2(isnan(TC2))=0;
        if Switch==1
            p=subplot(6,7,[i i+7]);
            p1=get(p,'position');
            if i==1;
                annotation('textbox',[p1(1)-p1(3)*0.3 p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],...
                    'string','a','EdgeColor','none','FontSize',18,'FontWeight','bold')
            end
            O=VecCorrCoef_v3_fast(TC1,TC2);
            imagescwithnan(center(Sq1.PRM.FrameBins),center(Sq1.PRM.FrameBins),(O'),jet,[1 1 1]);
            ax=gca;ax.XAxis.Color = CLR(i,:);ax.YAxis.Color = CLR(i,:);
            tempo=caxis;
            title([TEXTsmall{i} ' - Tuned'],'color',CLR(i,:))
            
        else
            p=subplot(6,7,[i+14 i+21]);
            p1=get(p,'position');
            if i==1;
                annotation('textbox',[p1(1)-p1(3)*0.3 p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],...
                    'string','b','EdgeColor','none','FontSize',18,'FontWeight','bold')
            end
            O=VecCorrCoef_v3_fast(TC1,TC2);
            imagescwithnan(center(Sq1.PRM.FrameBins),center(Sq1.PRM.FrameBins),(O'),jet,[1 1 1]);
            ax=gca;ax.XAxis.Color = [0.7 0.7 0.7];ax.YAxis.Color = [0.7 0.7 0.7];
            %caxis([tempo])
            title([TEXTsmall{i} ' - Untuned'],'color',[0.7 0.7 0.7])
            
        end
        yticks([400 800]);xticks([400 800])
        %colorbar('SouthOutside')
        %caxis([-1 1])
        %if i<4;    caxis([0.36 1]);    else;        caxis([0.6 1.05]);    end
        axis square
        if i==1
            line([1 600],[301 900],'color','k','linestyle','-.','linewidth',1);
            line([301 900],[1 600],'color','k','linestyle','-.','linewidth',1);
            if Switch==1
                cbr=colorbar('SouthOutside');
                ylabel(cbr,'Pop. overlap','FontSize',12);
            end
        end
        caxis([-0.4 0.8])
        %Type 2 The various diagonals
        clearvars Diags1 Diags2 Diags3
        for k=-899:899
            Diags1(k+900)=nanmean(diag(O,k));
            Diags2(k+900)=nanstd(diag(O,k));
            Diags3(k+900)=nanlength(diag(O,k));
        end
        %figure(2)
        p=subplot(6,7,[i+28 i+35]);
        p1=get(p,'position');
        if i==1;
            annotation('textbox',[p1(1)-p1(3)*0.3 p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],...
                'string','c','EdgeColor','none','FontSize',18,'FontWeight','bold')
        end
        if Switch==1
            shadedErrorBar(VecDiags,Diags1,Diags2./sqrt(Diags3),CLR(i,:));
            hold on
        else
            shadedErrorBar(VecDiags,Diags1,Diags2./sqrt(Diags3),[0.7 0.7 0.7]);
            hold on
        end
        axis tight;box off;ylim([-0.1 0.6])
        if i==1;
            xlabel('Diagonal number');ylabel('Avg. correlation');xticks([-450 -150 0 150 450])
            line([-300 -300],ylim,'color','k','linestyle','-.','linewidth',1);
            line([300 300],ylim,'color','k','linestyle','-.','linewidth',1);
        end
        
        if Switch==1
            DiagsOut(i,2,:)=Diags1;
        else
            DiagsOut_shf(i,2,:)=Diags1;
        end
        xlim([-450 450])
        %if i<4;ylim([0.5 1.25]);else;ylim([0.85 1.15]);end
        
    end
end

%% Compute the range of gaze
SList=(unique(DumpNS.ExpID(cell2mat(IDAll))))';
PupilPara=nan(length(D),3,2);
for ses=SList
    fileName=[D(ses).folder '\' D(ses).name];
    [ScP,~,~,ST,SI,TBins,RunFreeID,~,~,FrameID,EProp2]=...
        Allen_VVDot_Get_RunBehavior_Movies([fileName]...
        ,Sq1.PRM.TCutoff,Sq1.PRM.RVCut);
    ses
    PupilPara(ses,1,:)=[nanmin(EProp2([1:27000 27002:54000],1));nanmax(EProp2([1:27000 27002:54000],1))];
    PupilPara(ses,2,:)=[nanmin(EProp2([1:27000 27002:54000],3));nanmax(EProp2([1:27000 27002:54000],3))];
    PupilPara(ses,3,:)=[nanmin(EProp2([1:27000 27002:54000],4));nanmax(EProp2([1:27000 27002:54000],4))];
end

%% Supplement FR actual and that expected from eye movements

temp=A2([2 6 10 11 ]);
figure('Position',1e3*[ 0.2354    0.0330    1.2280    0.7472]);
%temp=temp([1 4]);
for i=1:4
    ses=DumpNS.ExpID(temp(i));
    UnitID=DumpNS.UnitID(temp(i));
    FileName=[D(ses).folder '\' D(ses).name];
    LocIN=[6,4,i;6,4,i+4;6,4,i+8;6,4,i+12;6,4,i+16;6,4,i+20;];
    if i==1;LabelBool=1;else; LabelBool=0;end
    FR_actC=Sq1.FR_act(temp(i),:)/delTCurrent;
    Allen_PlotMovieScatter_withEye_v5_MeanNorm(UnitID,FileName,FR_actC,Sq1.PRM.FrameBins,Sq1.PRM.TCutoff,Sq1.PRM.RVCut,60,LocIN,LabelBool);
end

%% Additional things about peak rate
PK_Cutoff=10;
clearvars ksp
for i=1:7
    idc=mintersect(Broad2,IDAll{i},Active2_sca,Active2);%find(Sc1.Zspa<2)
    tc1=Sq3.FR_act/Sq3.PRM.delTCurrent;
    tc2=Sc1.FR_act/Sq1.PRM.delTCurrent;
    clearvars pkRate mnRate
    for j=1:length(idc)
        pkRate(j,1)=nanmean(tc1(j,tc1(j,:)>prctile(tc1(j,:),100-PK_Cutoff)  ));
        pkRate(j,2)=nanmean(tc2(j,tc2(j,:)>prctile(tc2(j,:),100-PK_Cutoff)  ));
        
        mnRate(j,1)=nanmean(tc1(j,:));
        mnRate(j,2)=nanmean(tc2(j,:));
    end
    subplot(2,7,i)
    plot(pkRate(:,1),pkRate(:,2),'.','color',CLR(i,:))
    [ksp(i,1)]=signrank(log10(pkRate(:,1)),log10(pkRate(:,2)));
    [~,ksp(i,2)]=kstest2(log10(pkRate(:,1)),log10(pkRate(:,2)));
    xlim([0.5 100]);ylim([0.5 100]);
    set(gca,'xscale','log');set(gca,'yscale','log');
    box off;axis square;plot45Line([0.7 0.7 0.7])
    title(TEXTsmall{i},'color',CLR(i,:))
    if i==1;xlabel('Peak rate - Sequential (Hz)');ylabel('Scrambled (Hz)');end
    
    subplot(2,7,i+7)
    plot(mnRate(:,1),mnRate(:,2),'.','color',CLR(i,:))
    [ksp(i,3)]=signrank(log10(mnRate(:,1)),log10(mnRate(:,2)));
    [~,ksp(i,4)]=kstest2(log10(mnRate(:,1)),log10(mnRate(:,2)));
    xlim([0.5 100]);ylim([0.5 100]);
    set(gca,'xscale','log');set(gca,'yscale','log');
    box off;axis square;plot45Line([0.7 0.7 0.7])
    if i==1;xlabel('Mean rate - Sequential (Hz)');ylabel('Scrambled (Hz)');end
end

%_______________ Check the peak rate to scrambled movie with frames
%arranged as visually in order, or just by the time of presentation
TRWise_sca=load('Movie_FRStore_TrWise_AllTr_Sca.mat');
spar=nan(height(DumpNS),5);
idall=cell2mat(IDAll);idall=sortrows(idall);
for j=1:length(idall)
    i=idall(j);
    
    for tr=1:20
        tempA(tr,:)=TRWise_sca.TCDump_TrWise(i,tr,:);
        tempB(tr,:)=TRWise_sca.TCDump_TrWise(i,tr,MV.ShuffledFrameIDs);
        tempC(tr,:)=jmm_smooth_1d_cor(squeeze(TRWise_sca.TCDump_TrWise(i,tr,:)),Sq1.PRM.SMO,1);
        tempD(tr,:)=jmm_smooth_1d_cor(squeeze(TRWise_sca.TCDump_TrWise(i,tr,MV.ShuffledFrameIDs)),Sq1.PRM.SMO,1);
        for shifts=-10:10
            tempE(tr,shifts+11,:)=jmm_smooth_1d_cor...
                (squeeze(TRWise_sca.TCDump_TrWise(i,tr,circshift(MV.ShuffledFrameIDs,shifts))),Sq1.PRM.SMO,1);
        end
    end
    spar(i,1)=makeSparsity(nanmean(tempA));
    spar(i,2)=makeSparsity(nanmean(tempB));
    spar(i,3)=makeSparsity(nanmean(tempC));
    spar(i,4)=makeSparsity(nanmean(tempD));
    for shifts=1:21
        tempspar(shifts)=makeSparsity(squeeze(nanmean(tempE(:,shifts,:))));
    end
    spar(i,5)=max(tempspar);
    j
end

%% Subregions of V1 and CA1 plot
clc
load('Movie_Seq_Sca_CRCF_Rearranged_Feb2022_v2.mat');
[~,BestLatency]=nanmax(CRCF_dump_hr');
BestLatency(isnan(CRCF_dump_hr(:,1)))=nan;
load('Movie_Remapping_CRCF_Zscores_Feb2022.mat')

figure('Position',1e3*[0.2194    0.0794    1.2936    0.5768]);
LowCellCount=100;
LowCellCount_Tuned=50;

N_peaks=nan(height(DumpNS),1);
MedWidth=nan(height(DumpNS),1);
RatWidth=nan(height(DumpNS),1);
for i=1:height(DumpNS)
    N_peaks(i)=size(Sq5_hr.PkInfo{i},2);
    if N_peaks(i)>0
        MedWidth(i)=nanmedian(Sq5_hr.PkInfo{i}(3,:)+Sq5_hr.PkInfo{i}(4,:))*Sq1.PRM.delTCurrent;
        if N_peaks(i)>1
            RatWidth(i)=nanmax(Sq5_hr.PkInfo{i}(3,:)+Sq5_hr.PkInfo{i}(4,:))/...
                nanmin(Sq5_hr.PkInfo{i}(3,:)+Sq5_hr.PkInfo{i}(4,:));
        end
    end
end

%%%%%%%%%_____________Compute the tuning in different layers of V1
RegionIds={2:3;4:7};
LabelX={'Movie tuning (%)','Avg. number of peaks','Median peak width (sec)','Median peak width ratio','Scrambled tuning (%)','Scrambled tuned cell count','Remapping seq-sca'};
CLRRegion=[1 0 0;0 0 1];
for Region=1:2
    ID_Current=intersect(cell2mat(IDAll(RegionIds{Region})),Broad2);
    clearvars Name_subregion CountSubregions MovieTuning_Subregions
    
    ZCurrent=Sq1.Zspa(ID_Current);
    ZLater=Sc1.Zspa(ID_Current);
    ZRemapping=CRCF_z(ID_Current);
    Npeaks=N_peaks(ID_Current);
    MWidths=MedWidth(ID_Current);
    RWidths=RatWidth(ID_Current);
    a=DumpNS.RegionNumber(ID_Current);
    NumbersCurr=unique(a(~isnan(a)));
    for i=1:length(NumbersCurr)
        CountSubregions(i,1)=sum(a==NumbersCurr(i));
        MovieTuning_Subregions(i,1)=FracAbove2(ZCurrent(a==NumbersCurr(i)))*100;
        MovieTuning_Subregions(i,2)=nanmean(Npeaks(a==NumbersCurr(i) & ZCurrent>2));
        MovieTuning_Subregions(i,3)=nanmedian(MWidths(a==NumbersCurr(i) & ZCurrent>2));
        MovieTuning_Subregions(i,4)=nanmedian(RWidths(a==NumbersCurr(i) & ZCurrent>2));
        MovieTuning_Subregions(i,5)=FracAbove2(ZLater(a==NumbersCurr(i)))*100;
        MovieTuning_Subregions(i,6)=sum(a==NumbersCurr(i) & ZCurrent>2);
        MovieTuning_Subregions(i,7)=FracAbove2(ZRemapping(a==NumbersCurr(i)))*100;
        try
            Name_subregion{i,1}=SubregionTable.FullName{SubregionTable.ID==NumbersCurr(i)};
        catch
            Name_subregion{i,1}='Not_Identified';
        end
    end
    
    tempName=Name_subregion(CountSubregions>LowCellCount);
    for i=1:length(tempName)
        tempName{i}=strrep(tempName{i},'Primary visual area','V1');
        tempName{i}=strrep(tempName{i},'Posterior parietal association areas','PM');
        tempName{i}=strrep(tempName{i},'Anteromedial visual area','AM');
        tempName{i}=strrep(tempName{i},'Field','');
        tempName{i}=strrep(tempName{i},'Basic cell groups and regions','Unspecified');
        tempName{i}=strrep(tempName{i},'stratum','');
        %tempName{i}=strrep(tempName{i},'layer','');
        tempName{i}=strrep(tempName{i},'Subiculum','SUB');
        tempName{i}=strrep(tempName{i},'Dentate gyrus','DG');
    end
    for plt=1:7
        if Region==1;subplot_CSP(2,7,plt + (Region-1)*7,plt);else;subplot(2,7,plt + (Region-1)*7);end
        X = categorical(tempName);
        b=barh(X,MovieTuning_Subregions(CountSubregions>LowCellCount,plt));
        hold on
        if plt==1;line([5 5],ylim,'color','k');else;yticks([]);end
        if Region>-1;xlabel(LabelX{plt});end
        box off; axis square;axis tight
        b(1).FaceColor=CLRRegion(Region,:);b(1).EdgeColor=CLRRegion(Region,:);
    end
end

%% Supplement about the rearrangement stuff
load('Movie_Scra_Corrcoef_May1_2022.mat');
figure('Position',1e3*[0.4880    0.0802    1.0106    0.6818])

SmoLevel=Sq1.PRM.SMO*Sq5_hr.PRM.Upsample;
XBins=center(Sq5_hr_middleTrials.PRM.FrameBins);
ChosenID=[27555 19814 13315];%lgn 27555,am best -13315
LocBlock={[1 2 3 ],[5 6 7],[9 10 11]};
for i=1:3
    
    tc1=Sq5_hr_middleTrials.SpCnt_actual((ChosenID(i)),:)/SpCnt2FRScaler;
    tc2=Sc5_hr.SpCnt_actual((ChosenID(i)),:)/SpCnt2FRScaler;
    [bestCRCF,bestLatency]=max(CRCF_dump_hr(ChosenID(i),:));
    tc1_shifted=circshift(tc1,Timeshift(bestLatency));
    tc2_shifted=circshift(tc2,Timeshift(bestLatency));
    tc2_reassigned(ShuffledFrameIDs_hr)=tc2_shifted;
    p=subplot(4,4,LocBlock{i});
    
    plot(XBins,jmm_smooth_1d_cor_circ(tc2_shifted,SmoLevel,1),'color',[0.8 0.8 0.8],'linewidth',2);
    hold on;
    plot(XBins,jmm_smooth_1d_cor_circ(tc1_shifted,SmoLevel,1),'color',ColorSeq2,'linewidth',2);
    plot(XBins,jmm_smooth_1d_cor_circ(tc2_reassigned,SmoLevel,1),'color',ColorSca,'linewidth',2);
    box off
    axis tight
    xticks([400 800])
    if i==1;ylabel('Firing rate (Hz)');
        L=legend('Scrambled - as is','Sequential','Scrambled - rearranged');
        legend boxoff;L.ItemTokenSize(1)=8;
        p1=get(p,'position');
        annotation('textbox',[p1(1)-p1(3)/10 p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],'string','a',...
            'EdgeColor','none','FontSize',18,'FontWeight','bold')
    else
        if i==3
            xlabel('Frame #');
            
        end
    end
    title([TEXTsmall{i} ' r=' num2str(roundsd(bestCRCF,2)) ],'color',CLR(i,:));
    ax=gca;ax.YAxis.Color = CLR(i,:);ax.XAxis.Color = CLR(i,:);
end

p=subplot(4,4,[4 8]);
p1=get(p,'position');
annotation('textbox',[p1(1)-p1(3)/4 p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],'string','b',...
    'EdgeColor','none','FontSize',18,'FontWeight','bold')
for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq1.Zspa)),find(~isnan(Sc1.Zspa)),Active2 );
    [a,bb]=ecdf(CRCF_z(PYR));plot(bb,a,'linewidth',2,'color',CLR(i,:))
    hold on
end
box off;xlabel('Correlation(Episodic_{resp.},Scrambled_{resp.}) (z)')
xlim([-2 9]);ylabel('CDF');grid minor
line([2 2],ylim,'color','k','linestyle','--')
L=legend([TEXTsmall 'Significance threshold']); legend boxoff;L.ItemTokenSize(1)=8;

clearvars LowerFrac HigherFrac SummaryRearrProp DumpRearrProp
for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq1.Zspa)),Active2,TunedCells );
    ZCurrent=Sq1.Zspa(PYR);
    [bestCRCFValue_hr,LagBest_hr]=max(CRCF_dump_hr(PYR,:)');
    LagBest_hr_ms=Timeshift(LagBest_hr)*delTCurrent/10*1000;%10x upsample, and convert sec to ms by *1000
    RandJitter=(-0.5+rand(size(LagBest_hr))) * mean(diff(Timeshift))*delTCurrent/10*1000;
    LagBest_hr_ms=LagBest_hr_ms+ RandJitter;
    p=subplot(4,4,12)
    if i<4
        try
            [a,b]=ecdf(-LagBest_hr_ms(bestCRCFValue_hr>0.25));plot(b,a,'color',CLR(i,:),'linewidth',2)
            DumpRearrProp{i}=-LagBest_hr_ms(bestCRCFValue_hr>0.25);
            SummaryRearrProp(i,:)=[nanmean(DumpRearrProp{i}) nanmedian(DumpRearrProp{i}) nanstd(DumpRearrProp{i}) ...
                nanstd(DumpRearrProp{i})/sqrt(nanlength(DumpRearrProp{i}))];
        catch
            line([0 0],ylim,'linestyle','--','color',CLR(i,:),'linewidth',2)
        end
        hold on
    end
    DumpLags_seq_sca{i}=[-LagBest_hr_ms(bestCRCFValue_hr>0.25)];
    if i==3;
        L=legend(TEXTsmall(1:3)) ;legend boxoff;L.ItemTokenSize(1)=8; xlim([-10 150])
        box off;  title('r>0.25 cells only')
        xlabel('Latency for maximal correlation (ms)')
        ylabel('Probability');
        
        p1=get(p,'position');
        annotation('textbox',[p1(1)-p1(3)/2 p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],'string','d',...
            'EdgeColor','none','FontSize',18,'FontWeight','bold')
    end
    p=subplot(4,7,i+21);
    plot(-LagBest_hr_ms,bestCRCFValue_hr,'.','color',CLR(i,:),'markersize',4)
    axis tight;xlim([-350 350])
    if i==1
        xlabel('Latency for maximal correlation (ms)')
        ylabel('Max. correlation')
        
        p1=get(p,'position');
        annotation('textbox',[p1(1)-p1(3)/2 p1(2)+p1(4)*1.1 p1(3)/10 p1(4)/10],'string','c',...
            'EdgeColor','none','FontSize',18,'FontWeight','bold')
    end
    title(TEXTsmall{i},'color',CLR(i,:))
    ylim([0 0.6]);box off
    line(xlim,[0.25 0.25],'linestyle','--','color',[0.7 0.7 0.7])
    line([0 0],ylim,'linestyle','--','color',[0.7 0.7 0.7])
end

%% Simultaneuosly recorded cells
figure('Position',1e3*[    1.9602   -0.5126    1.0064    1.1840]);
CellIds=[ 38708 38911 38916 39329  ... DG Expt 44%
    26754 26760 26750 26758 ... CA3 ext 28
    13052 13243 13586 14401 ...  % CA1 Ses 22
    38120 38124 38252 38258; % SUB Expt 44
    ];
RegionCurrent=[4*ones(1,4) 5*ones(1,4) 6*ones(1,4) 7*ones(1,4)];
LabelBool=[1 zeros(1,15)];
ScatterSize=[ones(1,11) 2 ones(1,2) 1 2];
SMOCurrent=32;
TimeXAxisBool=0;
TimeXAxis=center(Sq5_hr.PRM.FrameBins)*Sq5_hr.PRM.delTCurrent;
for i=1:16
    idcurrent=CellIds(i);
    ses=DumpNS.ExpID(idcurrent);
    FileName=[D(ses).folder '\' D(ses).name];
    UnitID=DumpNS.UnitID(idcurrent);
    FR_actC=jmm_smooth_1d_cor_circ(Sq5_hr.SpCnt_actual(idcurrent,:)/SpCnt2FRScaler,SMOCurrent,1);
    PksAccepted=Sq5_hr.PkInfo{idcurrent};
    FR_at_peaks=FR_actC(nearestpoint(PksAccepted(2,:),center(Sq5_hr.PRM.FrameBins)));
    switch i
        case 1
            subplot_CSP(8,2,i,1);
        case 5
            subplot_CSP(8,2,i,2);
        case 9
            subplot_CSP(8,2,i,3);
        case 13
            subplot_CSP(8,2,i,4);
        otherwise
            subplot(8,2,i);
    end
    
    Movie_PlotScatter_psychedelic_v5_hr_Overlapping(UnitID,FileName,Sq1.PRM.TCutoff,Sq1.PRM.RVCut,LabelBool(i),CLR(RegionCurrent(i),:),ScatterSize(i),FR_actC,center(Sq5_hr.PRM.FrameBins));
    if rem(i,2)==0
        set(gca,'position',get(gca,'position')+[-0.04 0 0 0 ])
    end
end
%% Handful hippocampal cells with multiple peaks

figure('Position',1e3*[    1.9602   -0.5126    1.0064    0.840]);
CellIds=[ 12802  38708 43167 43282  ... DG Expt 44%
    53886  53889 54047 54052 ... CA3 ext 28
    23323 26730 38467 ...  % CA1 Ses 22
    25253; % SUB Expt 44
    ];
RegionCurrent=[4*ones(1,4) 5*ones(1,4) 6*ones(1,3) 7*ones(1,1)];
LabelBool=[1 zeros(1,15)];
ScatterSize=[ones(1,12)];
SMOCurrent=20;
TimeXAxisBool=0;
TimeXAxis=center(Sq5_hr.PRM.FrameBins)*Sq5_hr.PRM.delTCurrent;
for i=1:12
    idcurrent=CellIds(i);
    ses=DumpNS.ExpID(idcurrent);
    FileName=[D(ses).folder '\' D(ses).name];
    UnitID=DumpNS.UnitID(idcurrent);
    FR_actC=jmm_smooth_1d_cor_circ(Sq5_hr.SpCnt_actual(idcurrent,:)/SpCnt2FRScaler,SMOCurrent,1);
    switch i
        case 1
            subplot_CSP(6,2,i,1);
        otherwise
            subplot(6,2,i);
    end
    Movie_PlotScatter_psychedelic_v5_hr_Overlapping(UnitID,FileName,Sq1.PRM.TCutoff,Sq1.PRM.RVCut,LabelBool(i),CLR(RegionCurrent(i),:),ScatterSize(i),FR_actC,center(Sq5_hr.PRM.FrameBins));
    title(TEXTsmall{RegionCurrent(i)},'color',CLR(RegionCurrent(i),:));
    
end
%% Code for the circular shifts, simplified, because we dont care about the mean rates and whatnot
%{
Jumps=5;
Timeshift=-150:Jumps:150;
SMOLevelCurrent=Sq1.PRM.SMO*Sq5_hr.PRM.Upsample;
ShuffledFrameIDs_hr=nan(9000,1);
TopPerc=10;
for i=1:900
    currentID=MV.ShuffledFrameIDs(i);
    ShuffledFrameIDs_hr((i*10-9):(i*10))= currentID*10-4:currentID*10+5;
end
ShuffledFrameIDs_hr(ShuffledFrameIDs_hr==9001)=1;
ShuffledFrameIDs_hr(ShuffledFrameIDs_hr==9002)=2;
ShuffledFrameIDs_hr(ShuffledFrameIDs_hr==9003)=3;
ShuffledFrameIDs_hr(ShuffledFrameIDs_hr==9004)=4;
ShuffledFrameIDs_hr(ShuffledFrameIDs_hr==9005)=5;

BestLatency=nan(height(DumpNS),1);
CRCF_dump_hr=nan(height(DumpNS),length(Timeshift));
CRCF_dump_hr_asIs=nan(height(DumpNS),1);
CRCF_z=nan(height(DumpNS),1);
crcf_shf=nan(length(Timeshift),Sq1.PRM.nBoots);
CRCF_z_asis=nan(height(DumpNS),1);

clearvars RandomSequences
for n=1:Sq1.PRM.nBoots
    RandomSequences(:,n)=circshift(ShuffledFrameIDs_hr,randi(9000));
end

for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq1.Zspa)),Active2,TunedCells );
    
    for u=1:length(PYR)
        tc1=Sq5_hr_middleTrials.SpCnt_actual(PYR(u),:);
        tc1=jmm_smooth_1d_cor_circ(tc1,SMOLevelCurrent,1);
        tc2=Sc5_hr.SpCnt_actual(PYR(u),:);
        for tm=1:length(Timeshift)
            shiftCurrent=Timeshift(tm);
            tc1_shifted=circshift(tc1,(shiftCurrent));
            tc2_shifted=circshift(tc2,(shiftCurrent));
            
            for n=1:Sq1.PRM.nBoots
                tc2_reassignedRAND(RandomSequences(:,n))=tc2_shifted;
                tc2_reassignedRAND=jmm_smooth_1d_cor_circ(tc2_reassignedRAND,SMOLevelCurrent,1);
                crcf_shf(tm,n)=corrcoefwithNaNDiag(tc1_shifted,tc2_reassignedRAND);
            end
            
            tc2_reassigned(ShuffledFrameIDs_hr)=tc2_shifted;
            tc2_reassigned=jmm_smooth_1d_cor_circ(tc2_reassigned,SMOLevelCurrent,1);
            CRCF_dump_hr(PYR(u),tm)=corrcoefwithNaNDiag(tc1_shifted,tc2_reassigned);
            tc2_shifted=jmm_smooth_1d_cor_circ(tc2_shifted,SMOLevelCurrent,1);
        end
        
        act_best=max(CRCF_dump_hr(PYR(u),:));
        shf_best=max(crcf_shf);temp=zscore([shf_best act_best]);CRCF_z(PYR(u))=temp(end);
        
        u
    end
    i
end
save('Movie_Scra_Corrcoef_May1_2022.mat','CRCF_z_asis','CRCF_z','CRCF_dump_hr','BestLatency','SMOLevelCurrent','Timeshift','Jumps','ShuffledFrameIDs_hr','-v7.3')
%}

%% Need a function reporting the strain, age and gender of the mice used in the study
temp=unique(DumpNS.ExpID)';
k=1;
clearvars Age Gender Pheno Species
for ses=temp
    tempName=[D(ses).folder '\' D(ses).name];
    Age{k,1}=h5read(tempName, ['/general/subject/age']);
    Gender{k,1}=h5read(tempName, ['/general/subject/sex']);
    Pheno{k,1}=h5read(tempName, ['/general/subject/genotype']);
    Species{k,1}=h5read(tempName, ['/general/subject/species']);
    k=k+1;
end

%% Plot the correlation between Mega chronicity and isolation, rate instability and rate between 2 conditions
figure('position',1e3*[ 2.5890   -0.2142    1.1896    0.9968]);
clearvars tempCRF tempCRFP
load('Movie_PeakRelatedData_DumpPeakProps_May12_2022.mat');
load('TempMeanRates.mat');
StableInMovie=find(isbetween(MeanRates_runfree(:,1)./MeanRates_runfree(:,2),0.8,1.2));
StableInDrG=find(isbetween(MeanRates_runfree(:,3)./MeanRates_runfree(:,4),0.8,1.2) & isbetween(MeanRates_runfree(:,4)./MeanRates_runfree(:,5),0.8,1.2) & isbetween(MeanRates_runfree(:,3)./MeanRates_runfree(:,5),0.8,1.2) );
for i=1:7
    if i==1;subplot_CSP(6,7,i,1);else;subplot(6,7,i);end % Mean rate in movie vs blank
    idc=mintersect(Broad2,IDAll{i},StableInMovie);
    plot(0.5*(MeanRates_runfree(idc,1)+MeanRates_runfree(idc,2)),MeanRates_runfree(idc,6),'.','color',CLR(i,:),'markersize',4)
    xlim([0.005 100]);ylim([0.005 100]);xticks([0.1 10]);yticks([0.1 10])
    set(gca,'yscale','log');set(gca,'xscale','log');plot45Line([0.7 0.7 0.7 ]);grid on;box off
    %line([0.001 50],[0.005 250],'color','r');line([0.001 50],[0.0002 10],'color','r')
    if i==1;xlabel('FR_{movie} (Hz)')
        ylabel('FR_{blank} (Hz)')
    end
    axis square
    title(TEXTsmall{i},'color',CLR(i,:));
    [tempCRF(i,1),tempCRFP(i,1)]=corrcoefwithNaNDiag(log10(0.5*(MeanRates_runfree(idc,1)+MeanRates_runfree(idc,2))),log10(MeanRates_runfree(idc,6)));
    
    if i==1;subplot_CSP(6,7,i+7,2);else;subplot(6,7,i+7);end % FR vs the zscore of tuning
    idc=mintersect(Broad2,IDAll{i},Active2);
    plot(nanmean(Sq1.FR_shf_mean(idc,:),2)/delTCurrent,Sq1.Zspa(idc),'.','color',CLR(i,:),'markersize',4)
    xlim([0.5 100]);xticks([ 1 10 100]);ylim([-2 10.5]);yticks([-2 2 6 10])
    set(gca,'xscale','log');box off;
    if i==1;xlabel('Mean firing rate (Hz)')
        ylabel('Movie tuning (z)')
    end
    axis square
    %line([1 1],ylim,'linestyle','--','color',[0.7 0.7 0.7],'linewidth',2)
    [tempCRF(i,2),tempCRFP(i,2)]=corrcoefwithNaNDiag(log10(nanmean(Sq1.FR_shf_mean(idc,:),2)/delTCurrent),Sq1.Zspa(idc));
    
    if i==1;subplot_CSP(6,7,i+14,3);else;subplot(6,7,i+14);end % FR vs the number of peaks
    idc=mintersect(Broad2,IDAll{i},TunedCells,Active2);
    plot(nanmean(Sq1.FR_shf_mean(idc,:),2)/delTCurrent,DumpPeakProp{i,1}+ (2*rand(size(DumpPeakProp{i,1}))-1),'.','color',CLR(i,:),'markersize',4)
    xlim([0.5 100]);xticks([ 1 10 100]);if i<4;ylim([0 30]);yticks([0 10 20 30]);else;ylim([0 8]);yticks([0 4 8]);end
    set(gca,'xscale','log');box off;
    if i==1;xlabel('Mean firing rate (Hz)')
        ylabel('# Movie-fields')
    end
    axis square
    %line([1 1],ylim,'linestyle','--','color',[0.7 0.7 0.7],'linewidth',2)
    [tempCRF(i,3),tempCRFP(i,3)]=corrcoefwithNaNDiag(log10(nanmean(Sq1.FR_shf_mean(idc,:),2)/delTCurrent),DumpPeakProp{i,1}+ (2*rand(size(DumpPeakProp{i,1}))-1));
    
    if i==1;subplot_CSP(6,7,i+21,4);else;subplot(6,7,i+21);end % Megachronicity vs FR modulation in 2 blocks
    idc=mintersect(Broad2,IDAll{i},TunedCells,Active2);
    plot(nanmean(Sq1.FR_shf_mean(idc,:),2)/delTCurrent,DumpPeakProp{i,2},'.','color',CLR(i,:),'markersize',4)
    xlim([0.5 100]);xticks([ 1 10 100]);ylim([01 1000]);yticks([1 1000]);
    set(gca,'yscale','log');set(gca,'xscale','log');box off;
    if i==1;xlabel('Mean firing rate (Hz)')
        ylabel('Mega-scale index')
    end
    axis square
    %line([1 1],ylim,'linestyle','--','color',[0.7 0.7 0.7],'linewidth',2)
    [tempCRF(i,4),tempCRFP(i,4)]=corrcoefwithNaNDiag(log10(nanmean(Sq1.FR_shf_mean(idc,:),2)/delTCurrent),DumpPeakProp{i,2});
    
    if i==1;subplot_CSP(6,7,i+28,5);else;subplot(6,7,i+28);end % Megachronicity vs refractory violation
    idc=mintersect(Broad2,IDAll{i},TunedCells,Active2);
    plot(DumpNS.REF(idc),DumpPeakProp{i,2},'.','color',CLR(i,:),'markersize',4)
    ylim([01 1000]);xlim([0.001 100]);xticks([0.01 1 100]);yticks([1 1000])
    set(gca,'yscale','log');box off;set(gca,'xscale','log');%axis tight
    if i==1;ylabel('Mega-scale index')
        xlabel('Refractory violations index')
    end
    axis square
    [tempCRF(i,5),tempCRFP(i,5)]=corrcoefwithNaNDiag(DumpNS.REF(idc),DumpPeakProp{i,2});
    
    if i==1;subplot_CSP(6,7,i+35,6);else;subplot(6,7,i+35);end % Megachronicity vs isolation index
    idc=mintersect(Broad2,IDAll{i},TunedCells,Active2);
    plot(DumpNS.ISO(idc),DumpPeakProp{i,2},'.','color',CLR(i,:),'markersize',4)
    ylim([01 1000]);xlim([10 300]);xticks([1 10 100]);yticks([1 1000])
    set(gca,'yscale','log');set(gca,'xscale','log');box off;
    if i==1;ylabel('Mega-scale index')
        xlabel('Isolation index')
    end
    axis square
    [tempCRF(i,6),tempCRFP(i,6)]=corrcoefwithNaNDiag(DumpNS.ISO(idc),DumpPeakProp{i,2});
end

%% Create an Extended Data figure where I remove the SWRs - copy the format of RunFree figure creation
load('Allen_SWRandThetaSupport_v1.mat');
SMOVariable=[1 1 1 2 1 2 2;1 1 1 2 2 2 2;];
figure('position',[510.6000  -22.2000  989.6000  775.2000]);
SPLoc=[1:7];
UseAllData=0;
ReplaceSWRByRun=1;
A2Current=A2;
A2Current(4,1)=27048;
for j=1:size(A2Current,1)
    for i=1:1%size(A2,2)
        idcurrent=A2Current(j,i);
        ses=DumpNS.ExpID(idcurrent);
        FileName=[D(ses).folder '\' D(ses).name];
        UnitID=DumpNS.UnitID(idcurrent);
        FR_actC=Sq6.FR_act(idcurrent,:)/delTCurrent;
        FR_actC=jmm_smooth_1d_cor_circ(FR_actC,Sq1.PRM.SMO*SMOVariable(i,j),1);
        %FR_actC=jmm_smooth_1d_cor_circ(FR_actC,Sq1.PRM.SMO,1);
        s1=subplot(3,3,SPLoc(j));
        LocIN=[s1.Position];if i==1; LocIN= LocIN+[0 0.0 0 0];end
        if  j==1
            LabelBool=1;
            annotation('textbox',[LocIN(1)-LocIN(3)*0.20 LocIN(2)+LocIN(4)*1.11 LocIN(3)/10 LocIN(4)/10],...
                'string','a','EdgeColor','none','FontSize',18,'FontWeight','bold')
        else
            LabelBool=0;
        end
        l1= Movie_PlotScatter_PSTH_psychedelic_v4_noStd(UnitID,FileName,Sq6.PRM.TCutoff,Sq6.PRM.RVCut,Sq6.PRM.FrameBins,...
            FR_actC,LocIN,LabelBool,CLR(j,:),UseAllData,ReplaceSWRByRun,DumpSesWise{ses},Sq6.PRM.SWRThreshold,Sq6.PRM.TCutoff_SWR);
        set(gcf,'CurrentAxes',l1)
        title(TEXTsmall{j},'color',CLR(j,:));
        %title([num2str(idcurrent) '-' DumpNS.Region(idcurrent)]);
    end
end

Zbins=linspace(-2.5,8,11);
clearvars store
s1=subplot(3,3,8);
for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq6.Zspa)),Active2 );store(i,1:3)=[nanlength(PYR) nansum(Sq6.Zspa(PYR)>2) FracAbove2(Sq6.Zspa(PYR))] ;
    ALL=(IDAll{i} );
    [a,b]=ecdf(Sq6.Zspa(PYR),'function','survivor');
    plot(b,a,'color',CLR(i,:),'linewidth',2)
    hold on
end
PYR=mintersect(Broad2,cell2mat(IDAll),find(~isnan(Sq2.Zspa)),Active2 );
[a,b]=ecdf(Sq6.Zspa_shf(PYR),'function','survivor');
plot(b,a,'color','k','linewidth',0.5,'linestyle','--')
hold on
xlim([-2 10]);xlabel('Movie tuning (z)')
line([2 2],ylim,'color',[0 0 0]);
legend(TEXTsmall);legend boxoff;L=legend;L.ItemTokenSize(1)=8;
ylabel('Probability')
box off
grid on
LocIN=[s1.Position];
annotation('textbox',[LocIN(1)-LocIN(3)*0.20 LocIN(2)+LocIN(4)*1.11 LocIN(3)/10 LocIN(4)/10],...
    'string','b','EdgeColor','none','FontSize',18,'FontWeight','bold')

clearvars storeAll storeSta storeEq ks_pvalue
for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq1.Zspa)),Active2 );
    storeAll(i,1:3)=[nanlength(PYR) nansum(Sq1.Zspa(PYR)>2) FracAbove2(Sq1.Zspa(PYR))] ;
    PYR2d=mintersect(Broad2,IDAll{i},find(~isnan(Sq6.Zspa)),Active2 );
    storeSta(i,1:3)=[nanlength(PYR2d) nansum(Sq6.Zspa(PYR2d)>2) FracAbove2(Sq6.Zspa(PYR2d))] ;
    PYR2d=mintersect(Broad2,IDAll{i},find(~isnan(Sq6_dash.Zspa)),Active2 );
    storeEq(i,1:3)=[nanlength(PYR2d) nansum(Sq6_dash.Zspa(PYR2d)>2) FracAbove2(Sq6_dash.Zspa(PYR2d))] ;
    [~,ks_pvalue(i)]=kstest2(Sq6.Zspa(PYR),Sq6_dash.Zspa(PYR2d));
end
s1=subplot(3,3,9);
b=bar([storeAll(:,end) storeSta(:,end) storeEq(:,end)]*100);
box off;
for i=1:3
    b(i).EdgeColor=b(i).FaceColor;
end
box off
ylabel('Movie tuning (%)')
xticks([1:7])
xticklabels(TEXTsmall)
xtickangle(20)
line(xlim,[2.25 2.25],'color',[0 0 0],'linestyle','--')
legend('All data','SWR removed data','Equivalent subsample','Chance level')
L=legend;legend boxoff;L.ItemTokenSize(1)=8;
LocIN=[s1.Position];
annotation('textbox',[LocIN(1)-LocIN(3)*0.20 LocIN(2)+LocIN(4)*1.11 LocIN(3)/10 LocIN(4)/10],...
    'string','c','EdgeColor','none','FontSize',18,'FontWeight','bold')
%% Create the additional plot about Pupil area scatter with color depicting the area, or similar for theta later
load('Allen_SWRandThetaSupport_v1.mat');
SMOVariable=[1 1 1 2 1 2 2;1 1 1 2 2 2 2;];
figure('position',1e3*[0.7114    1.1138    1.2488    0.5992]);
clearvars FR_actC
SPLoc=[1:7];
UseAllData=1;
A2Current=A2;
A2Current(4,1)=27048;
A2Current_last=[A2Current(1,1) A2Current(2,2) A2Current(3,1) A2Current(4,2) A2Current(5,1) A2Current(6,1) A2Current(7,2)];
A2Current_last=A2Current_last([2 6]);
A2Current_last=A2Current_last';
TEXTCurrent=TEXTsmall([2 6]);
CLRCurrent=CLR([2 6],:);
for j=1:size(A2Current_last,1)
    for i=1:1%size(A2,2)
        idcurrent=A2Current_last(j,i);
        ses=DumpNS.ExpID(idcurrent);
        FileName=[D(ses).folder '\' D(ses).name];
        UnitID=DumpNS.UnitID(idcurrent);
        %________________ Setup for later
        tempName=[D(ses).folder '\' D(ses).name];
        [ScP,TStart,TEnd,ST,SI,TBins,RunFreeID,RunForW,RV,FrameID,EProp2]=...
            Allen_VVDot_Get_RunBehavior_Movies(tempName...
            ,Sq1.PRM.TCutoff,Sq1.PRM.RVCut);
        delTCurrent=nanmedian(diff(TBins));
        RepCount=ceil((1:54000)/900);
        MovieDataID=zeros(length(center(TBins)),1);MovieDataID([1:27000 27002:length(RunFreeID)])=1;
        RunFreeID=RunFreeID(MovieDataID==1);
        
        %________________ Pupil area first
        temp=EProp2(:,1);
        PA=temp(MovieDataID==1);
        FR_actC(:,1)=Sq7.FR_act(idcurrent,:)/delTCurrent;
        FR_actC(:,1)=jmm_smooth_1d_cor_circ(FR_actC(:,1),Sq7.PRM.SMO*SMOVariable(i,j),1);
        FR_actC(:,2)=Sq7_dash.FR_act(idcurrent,:)/delTCurrent;
        FR_actC(:,2)=jmm_smooth_1d_cor_circ(FR_actC(:,2),Sq7_dash.PRM.SMO*SMOVariable(i,j),1);
        s1=subplot(4,3,j+6);%s1=subplot(3,7,j);
        LocIN=[s1.Position];if i==1; LocIN= LocIN+[0 0.0 0 0];end
        if  j==1
            LabelBool=1;
            annotation('textbox',[LocIN(1)-LocIN(3)*0.20 LocIN(2)+LocIN(4)*1.11 LocIN(3)/10 LocIN(4)/10],...
                'string','c','EdgeColor','none','FontSize',18,'FontWeight','bold')
        else
            LabelBool=0;
        end
        l1= Movie_PlotScatter_PSTH_psychedelic_v5_ColorDots(UnitID,FileName,Sq7.PRM.TCutoff,Sq7.PRM.RVCut,Sq7.PRM.FrameBins,...
            FR_actC,LocIN,LabelBool,CLRCurrent(j,:),UseAllData,PA,0);
        set(gcf,'CurrentAxes',l1)
        title(TEXTCurrent{j},'color',CLRCurrent(j,:));
        
        %________________ Theta power next
        ThetaPowerFull=(DumpSesWise{ses}([4],:))';
        TimeForTheta=DumpSesWise{ses}(1,:);
        temp=interp1(TimeForTheta,ThetaPowerFull,center(TBins));
        TP=temp(MovieDataID==1);
        FR_actC(:,1)=Sq8.FR_act(idcurrent,:)/delTCurrent;
        FR_actC(:,1)=jmm_smooth_1d_cor_circ(FR_actC(:,1),Sq8.PRM.SMO*SMOVariable(i,j),1);
        FR_actC(:,2)=Sq8_dash.FR_act(idcurrent,:)/delTCurrent;
        FR_actC(:,2)=jmm_smooth_1d_cor_circ(FR_actC(:,2),Sq8_dash.PRM.SMO*SMOVariable(i,j),1);
        s1=subplot(4,3,j+9);%s1=subplot(3,7,j+7);
        LocIN=[s1.Position];if i==1; LocIN= LocIN+[0 0.0 0 0];end
        if  j==1
            LabelBool=1;
            annotation('textbox',[LocIN(1)-LocIN(3)*0.20 LocIN(2)+LocIN(4)*1.11 LocIN(3)/10 LocIN(4)/10],...
                'string','d','EdgeColor','none','FontSize',18,'FontWeight','bold')
        else
            LabelBool=0;
        end
        [l1,l2,Clrmp]= Movie_PlotScatter_PSTH_psychedelic_v5_ColorDots(UnitID,FileName,Sq7.PRM.TCutoff,Sq7.PRM.RVCut,Sq7.PRM.FrameBins,...
            FR_actC,LocIN,LabelBool,CLRCurrent(j,:),UseAllData,TP,0);
        set(gcf,'CurrentAxes',l1)
    end
end

clearvars storePA
for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq7.Zspa) & ~isnan(Sq7_dash.Zspa)),Active2 );
    storePA(i,1:3)=[nanlength(PYR) nansum(Sq7.Zspa(PYR)>2) FracAbove2(Sq7.Zspa(PYR))] ;
    storePA(i,4:6)=[nanlength(PYR) nansum(Sq7_dash.Zspa(PYR)>2) FracAbove2(Sq7_dash.Zspa(PYR))] ;
    storePA(i,7:8)=[FracAbove2(Sq7_dash.Zspa_shf(PYR)) FracAbove2(Sq7_dash.Zspa_shf(PYR))];
    [~,ks_pvalue(1,i)]=kstest2(Sq7.Zspa(PYR),Sq7_dash.Zspa(PYR));
    [~,ks_pvalue(2,i)]=kstest2(Sq7.Zspa(PYR),Sq7.Zspa_shf(PYR));
    [~,ks_pvalue(3,i)]=kstest2(Sq7_dash.Zspa(PYR),Sq7_dash.Zspa_shf(PYR));
end
s1=subplot(4,3,9);
b=bar([storePA(:,[3 6:8])]*100);
box off;
b(1).FaceColor=Clrmp(1,:);b(3).FaceColor=[0.7 0.7 0.7];
b(2).FaceColor=Clrmp(end,:);b(4).FaceColor=[0.7 0.7 0.7];

for i=1:2
    b(i).EdgeColor=b(i).FaceColor;
    b(i+2).EdgeColor=b(i).FaceColor;
end
box off
ylabel('Movie tuning (%)')
xticks([1:7])
xticklabels(TEXTsmall)
xtickangle(20)
line(xlim,[2.25 2.25],'color',[0 0 0],'linestyle','--')
legend('Pupil constriction (bottom 50 percentile)','Pupil dilation (top 50 percentile)','Chance level','Chance level')
L=legend;legend boxoff;L.ItemTokenSize(1)=8;
LocIN=[s1.Position];
annotation('textbox',[LocIN(1)-LocIN(3)*0.20 LocIN(2)+LocIN(4)*1.11 LocIN(3)/10 LocIN(4)/10],...
    'string','e','EdgeColor','none','FontSize',18,'FontWeight','bold')
%______________________ Repeat for theta

clearvars storeTP
for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq8.Zspa) & ~isnan(Sq8_dash.Zspa)),Active2 );
    storeTP(i,1:3)=[nanlength(PYR) nansum(Sq8.Zspa(PYR)>2) FracAbove2(Sq8.Zspa(PYR))] ;
    storeTP(i,4:6)=[nanlength(PYR) nansum(Sq8_dash.Zspa(PYR)>2) FracAbove2(Sq8_dash.Zspa(PYR))] ;
    storeTP(i,7:8)=[FracAbove2(Sq8_dash.Zspa_shf(PYR)) FracAbove2(Sq8_dash.Zspa_shf(PYR))];
   [~,ks_pvalue(1,i)]=kstest2(Sq8.Zspa(PYR),Sq8_dash.Zspa(PYR));
    [~,ks_pvalue(2,i)]=kstest2(Sq8.Zspa(PYR),Sq8.Zspa_shf(PYR));
    [~,ks_pvalue(3,i)]=kstest2(Sq8_dash.Zspa(PYR),Sq8_dash.Zspa_shf(PYR));
end
s1=subplot(4,3,12);
b=bar([storeTP(:,[3 6:8])]*100);
box off;
b(1).FaceColor=Clrmp(1,:);b(3).FaceColor=[0.7 0.7 0.7];
b(2).FaceColor=Clrmp(end,:);b(4).FaceColor=[0.7 0.7 0.7];

for i=1:2
    b(i).EdgeColor=b(i).FaceColor;
    b(i+2).EdgeColor=b(i).FaceColor;
end
box off
ylabel('Movie tuning (%)')
xticks([1:7])
xticklabels(TEXTsmall)
xtickangle(20)
line(xlim,[2.25 2.25],'color',[0 0 0],'linestyle','--')
legend('Low theta power (bottom 50 percentile)','High theta power (top 50 percentile)','Chance level','Chance level')
L=legend;legend boxoff;L.ItemTokenSize(1)=8;
LocIN=[s1.Position];
annotation('textbox',[LocIN(1)-LocIN(3)*0.20 LocIN(2)+LocIN(4)*1.11 LocIN(3)/10 LocIN(4)/10],...
    'string','f','EdgeColor','none','FontSize',18,'FontWeight','bold')

%__________________ Plotting eg cells at the end
%CellIdsRunStaCompare=[49294;49857;49982;26754;26760;28273;];
%CellIdsRunStaCompare=[14808;14411;15075;26754;26760;28273;];
CellIdsRunStaCompare=[14808;14411;38138;26754;51882;28273];

RegionCurrent=[4;6;7;5;6;7;];
LabelBool=[1 zeros(1,5)];
ScatterSize=[ones(1,6)];
SMOVariable=[5;5;2;2 ;5 ;2;];
for i=1:6
    idcurrent=CellIdsRunStaCompare(i);
    ses=DumpNS.ExpID(idcurrent);
    FileName=[D(ses).folder '\' D(ses).name];
    UnitID=DumpNS.UnitID(idcurrent);
        FR_actC=Sq1.FR_act(idcurrent,:)/delTCurrent;
        FR_actC=jmm_smooth_1d_cor_circ(FR_actC,Sq1.PRM.SMO*SMOVariable(i),1);
    switch i
        case 1
            [~,s1]=subplot_CSP(4,3,i,1);
        case 4
            [~,s1]=subplot_CSP(4,3,i,2);
        otherwise
            s1=subplot(4,3,i);
            s1=s1.Position;
    end
    LocIN=[s1];
    [~,~,~,ST,SI,TBins,RunFreeID,~,RV,FrameID]=Allen_VVDot_Get_RunBehavior_Movies([FileName],Sq1.PRM.TCutoff,Sq1.PRM.RVCut);
    %Movie_PlotScatter_PSTH_psychedelic_v8_RVColoredSpikes(UnitID,FileName,Sq1.PRM.TCutoff,Sq1.PRM.RVCut,(Sq1.PRM.FrameBins),FR_actC,LabelBool,CLR(RegionCurrent(i),:),1,0)
       %l1= Movie_PlotScatter_PSTH_Seq_v4_Overlapping(UnitID,FileName,Sq1.PRM.TCutoff,Sq1.PRM.RVCut,Sq1.PRM.FrameBins,FR_actC,0,0,LocIN,LabelBool,CLR(RegionCurrent(i),:),1);
    
    Movie_PlotScatter_psychedelic_v5_hr_Overlapping(UnitID,FileName,Sq1.PRM.TCutoff,Sq1.PRM.RVCut,LabelBool(i),CLR(RegionCurrent(i),:),ScatterSize(i),FR_actC,center(Sq1.PRM.FrameBins));
        
    title([TEXTsmall{RegionCurrent(i)} ', %_{stationary time}= ' num2str(round(100*sum(RunFreeID)/length(RunFreeID)))],'color',CLR(RegionCurrent(i),:));
end
%% Another extended data figure, comparing depth of modulation and sparsity etc
Zmetrics=load('MovieDifferentTuningMetrics_Jan28_2023.mat');
figure('position',1e3*[2.1954    0.1010    1.3736    0.5432]);
clearvars storeMetrics KSTest CRCF
for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq1.Zspa)),Active2 );
    %{
    if i==1;subplot_CSP(7,5,(i-1)*5+1,1);else;subplot(7,5,(i-1)*5+1);end
    plot(nanmean(Sq1.FR_act(PYR,:)'),Zmetrics.SparOut(PYR,2),'.','color',ColorSca,'markersize',1);ylim([0 0.8])
    set(gca,'xscale','log');
    yyaxis right
    plot(nanmean(Sq1.FR_act(PYR,:)'),Zmetrics.SparOut(PYR,4),'.','color',ColorSeq,'markersize',1);ylim([-3 11])
    
    if i==1;subplot_CSP(3,7,i,1);else;subplot(3,7,i);end
    plot(Zmetrics.SparOut(PYR,4),Zmetrics.DepthOfMod(PYR,4),'.','color',CLR(i,:),'markersize',2);xlim([-3 11]);ylim([-3 11]);axis square;plot45Line([0.7 0.7 0.7]);box off
    line([2 2],ylim,'linestyle','--','color',[0.7 0.7 0.7]);line(xlim,[2 2],'linestyle','--','color',[0.7 0.7 0.7])
    if i==1;xlabel('Sparsity (z)');ylabel('Depth of modulation (z)');end
    title(TEXTsmall{i},'color',CLR(i,:))
    
    if i==1;subplot_CSP(3,7,i+7,2);else;subplot(3,7,i+7);end
    plot(Zmetrics.SparOut(PYR,4),Zmetrics.MutualInfo(PYR,4),'.','color',CLR(i,:),'markersize',2);xlim([-3 11]);ylim([-3 11]);axis square;plot45Line([0.7 0.7 0.7]);box off
    line([2 2],ylim,'linestyle','--','color',[0.7 0.7 0.7]);line(xlim,[2 2],'linestyle','--','color',[0.7 0.7 0.7])
    if i==1;xlabel('Sparsity (z)');ylabel('Mutual information (z)');end
    %}
    storeMetrics(i,:)=[FracAbove2(Zmetrics.SparOut(PYR,4)) FracAbove2(Zmetrics.DepthOfMod(PYR,4)) FracAbove2(Zmetrics.MutualInfo(PYR,4)) ...
        FracAbove2(Zmetrics.SparOut(PYR,3)) FracAbove2(Zmetrics.DepthOfMod(PYR,3)) FracAbove2(Zmetrics.MutualInfo(PYR,3))];
    
    [~,KSTest(i,1)]=kstest2(Zmetrics.SparOut(PYR,4),Zmetrics.DepthOfMod(PYR,4));[CRCF(i,1),CRCF(i,2)]=corrcoefwithNaNDiag(Zmetrics.SparOut(PYR,4),Zmetrics.DepthOfMod(PYR,4));
    [~,KSTest(i,2)]=kstest2(Zmetrics.SparOut(PYR,4),Zmetrics.MutualInfo(PYR,4));[CRCF(i,3),CRCF(i,4)]=corrcoefwithNaNDiag(Zmetrics.SparOut(PYR,4),Zmetrics.MutualInfo(PYR,4));
    
    [~,KSTest(i,3)]=kstest2(Zmetrics.DepthOfMod(PYR,4),Zmetrics.DepthOfMod(PYR,3));
    [~,KSTest(i,4)]=kstest2(Zmetrics.MutualInfo(PYR,4),Zmetrics.MutualInfo(PYR,3));
    [~,KSTest(i,5)]=kstest2(Zmetrics.SparOut(PYR,4),Zmetrics.SparOut(PYR,3));
end

%subplot_CSP(3,2,5,3);
b=bar([storeMetrics]*100);
box off;
b(1).FaceColor=ColorSeq2;b(4).FaceColor=[0.7 0.7 0.7];
b(2).FaceColor=ColorSca;b(5).FaceColor=[0.7 0.7 0.7];
b(3).FaceColor=co('ucla gold');b(6).FaceColor=[0.7 0.7 0.7];

for i=1:3
    b(i).EdgeColor=b(i).FaceColor;
    b(i+2).EdgeColor=b(i).FaceColor;
end
box off
ylabel('Movie tuning (%)')
xticks([1:7])
xticklabels(TEXTsmall)
xtickangle(20)
line(xlim,[2.25 2.25],'color',[0 0 0],'linestyle','--')
legend('z-scored sparsity','z-scored depth of modulation','z-scored mutual information','Chance level','Chance level','Chance level')
L=legend;legend boxoff;L.ItemTokenSize(1)=8;

%% Compare the zscores in actual movie and the rearranged one
figure('position',1e3*[0.5674    1.1746    1.3736    0.5432]);
ColorSca2=[0.2 0.2 0.6];
Zscore_SeqScaRearr=load('Movie_Zcompare_SeqScaRearr_Jan30_2023.mat');
clearvars storeMetrics KSTest CRCF
for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq1.Zspa)),Active2 );
    
    [CRCF(i,1),CRCF(i,2)]=corrcoefwithNaNDiag(Zscore_SeqScaRearr.SparOut(PYR,1),Zscore_SeqScaRearr.SparOut(PYR,5));
    [~,KSTest(i,1)]=kstest2(Zscore_SeqScaRearr.SparOut(PYR,1),Zscore_SeqScaRearr.SparOut(PYR,5));
    
    [CRCF(i,3),CRCF(i,4)]=corrcoefwithNaNDiag(Zscore_SeqScaRearr.SparOut(PYR,1),Zscore_SeqScaRearr.SparOut(PYR,9));
    [~,KSTest(i,2)]=kstest2(Zscore_SeqScaRearr.SparOut(PYR,1),Zscore_SeqScaRearr.SparOut(PYR,9));
    [~,KSTest_withShf(i,1)]=kstest2(Zscore_SeqScaRearr.SparOut(PYR,10),Zscore_SeqScaRearr.SparOut(PYR,9));
    
    [CRCF(i,5),CRCF(i,6)]=corrcoefwithNaNDiag(Zscore_SeqScaRearr.SparOut(PYR,5),Zscore_SeqScaRearr.SparOut(PYR,9));
    [~,KSTest(i,3)]=kstest2(Zscore_SeqScaRearr.SparOut(PYR,5),Zscore_SeqScaRearr.SparOut(PYR,9));
    
    storeMetrics(i,1:4)=[FracAbove2(Zscore_SeqScaRearr.SparOut(PYR,1)) FracAbove2(Zscore_SeqScaRearr.SparOut(PYR,5)) FracAbove2(Zscore_SeqScaRearr.SparOut(PYR,9)) FracAbove2(Zscore_SeqScaRearr.SparOut(PYR,10))];
end

b=bar(1:3,[storeMetrics(1:3,:)]*100);
box off;
b(1).FaceColor=ColorSeq2;b(1).EdgeColor=ColorSeq2;
b(2).FaceColor=ColorSca;b(2).EdgeColor=ColorSca;
b(3).FaceColor=ColorSca2;b(3).EdgeColor=ColorSca2;
b(4).FaceColor=[0.7 0.7 0.7];b(4).EdgeColor=[0.7 0.7 0.7];
ylabel('Visual tuning(%)')
line([1 3],[2.25 2.25],'color',[0 0 0],'linestyle','--')
ylim([0 100])
yyaxis right

b=bar(4:7,[storeMetrics(4:7,:)]*100);
line([4 7],[2.25 2.25],'color',[0 0 0],'linestyle','--')
b(1).FaceColor=ColorSeq2;b(1).EdgeColor=ColorSeq2;
b(2).FaceColor=ColorSca;b(2).EdgeColor=ColorSca;
b(3).FaceColor=ColorSca2;b(3).EdgeColor=ColorSca2;
b(4).FaceColor=[0.7 0.7 0.7];b(4).EdgeColor=[0.7 0.7 0.7];
xticks([1:7])
xticklabels([TEXTsmall])
xtickangle(20)
ax=gca;ax.YAxis(1).Color ='k';ax.YAxis(2).Color ='k';
line([3.5 3.5],ylim,'linestyle','-.')
box off
ylabel('Hippocampal tuning (%)')
legend('Continuous (middle 20 trials)','Scrambled (as is)','Scrambled (rearranged)','Chance level')
ylim([0 35])

%% For each cell, we want to consider the p value of the t test between responses at frame 1 vs frame 900, and a comparable pair of frames with lots of visual change, somewhere away
Pval_ofT=nan(height(DumpNS),6);
for i=1:7
    PYR=mintersect(Broad2,IDAll{i},find(~isnan(Sq1.Zspa)),Active2 );
   for j=1:length(PYR)
       frs1=[squeeze(TRW_seq.TCDump_TrWise(PYR(j),:,900)) ;squeeze(TRW_seq.TCDump_TrWise(PYR(j),:,1))];[~,Pval_ofT(PYR(j),1)]=ttest2(frs1(1,:),frs1(2,:));
       frs1=[squeeze(TRW_seq.TCDump_TrWise(PYR(j),:,1)) ;squeeze(TRW_seq.TCDump_TrWise(PYR(j),:,2))];[~,Pval_ofT(PYR(j),2)]=ttest2(frs1(1,:),frs1(2,:));
       frs1=[squeeze(TRW_seq.TCDump_TrWise(PYR(j),:,2)) ;squeeze(TRW_seq.TCDump_TrWise(PYR(j),:,3))];[~,Pval_ofT(PYR(j),3)]=ttest2(frs1(1,:),frs1(2,:));
       
       frs1=[squeeze(TRW_seq.TCDump_TrWise(PYR(j),:,155)) ;squeeze(TRW_seq.TCDump_TrWise(PYR(j),:,156))];[~,Pval_ofT(PYR(j),4)]=ttest2(frs1(1,:),frs1(2,:));
       frs1=[squeeze(TRW_seq.TCDump_TrWise(PYR(j),:,156)) ;squeeze(TRW_seq.TCDump_TrWise(PYR(j),:,157))];[~,Pval_ofT(PYR(j),5)]=ttest2(frs1(1,:),frs1(2,:));
       frs1=[squeeze(TRW_seq.TCDump_TrWise(PYR(j),:,157)) ;squeeze(TRW_seq.TCDump_TrWise(PYR(j),:,158))];[~,Pval_ofT(PYR(j),6)]=ttest2(frs1(1,:),frs1(2,:));
       [i,j]
   end
   subplot(2,4,i)
   [a,b]=ecdf(Pval_ofT(PYR,3));plot(b,a,'r','linewidth',2)
   hold on
   [a,b]=ecdf(Pval_ofT(PYR,6));plot(b,a,'b','linewidth',2)
   line([0.05 0.05],ylim,'color','k','linewidth',2,'linestyle','--')
   title(TEXTsmall{i})
   if i==1;xlabel('Pvalue of the ttest between neural responses from 60 trials');
       legend('1st vs 900th frame',' Another instance of low F2F correlation in the middle of the movie','0.05 line');end
end