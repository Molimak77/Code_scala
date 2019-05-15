## ================================================= ##
#                 Le projet Terry-Fox                ##
## ================================================= ##
rm(list=ls())
library(nloptr)
ch1<-c("~/Google Drive/VincentMoi/1ProjEV/Proj_CPCBN_Terry_fox/Dos_data/Data_RData/datMarq.csv")
ch2<-c("~/Google Drive/VincentMoi/1ProjEV/Proj_CPCBN_Terry_fox/Dos_data/Data_RData/datclinic.csv")

datmarq1<-read.csv(ch1,header = TRUE)
datclini<-read.csv(ch2,header = TRUE)
datwork<-merge(datmarq1,datclini,by="NID")

#datmarq<-datwork[datwork$datBRC!=0,]

datmarq<-datwork
# data management des données



datmarq$BCRMonthsorLastContact<-NULL
datmarq$MetsSite<-NULL




datmarq$Radiation01<-NULL
datmarq$Hormonotherapy01<-NULL
datmarq$Chemo01<-NULL
datmarq$MonthstoBoneMets<-NULL
datmarq$MonthstoBCR<-NULL

datmarq$BRC1<-factor(datmarq$BRC,labels=c("No","Yes"))
datmarq$StatusCRPC01<-factor(datmarq$StatusCRPC01,labels=c("No","Yes"))
datmarq$DeathOverall01<-factor(datmarq$DeathOverall01,labels=c("No","Yes"))
datmarq$PSA<-as.numeric(datmarq$PSA)
datmarq$pTNM1234<-factor(datmarq$pTNM1234)
datmarq$Margin01<-factor(ifelse(datmarq$Margin01==".",NA,datmarq$Margin01))

# reserve base

datmarqR<-datmarq

# fonctions utiles

source('~/Google Drive/VincentMoi/1ProjEV/StatModeling.v.01/abr_descft.R', chdir = TRUE)
source('~/Google Drive/VincentMoi/1ProjEV/StatModeling.v.01/Desc_fct.R', chdir = TRUE)
source('~/Google Drive/VincentMoi/1ProjEV/StatModeling.v.01/Desc_Qte.R', chdir = TRUE)
source('~/Google Drive/VincentMoi/1ProjEV/StatModeling.v.01/Dsc_couv.R', chdir = TRUE)
source('~/Google Drive/VincentMoi/1ProjEV/StatModeling.v.01/fct_NA.R', chdir = TRUE)
source('~/Google Drive/VincentMoi/1ProjEV/StatModeling.v.01/Fct_Rscox.R', chdir = TRUE)
source('~/Google Drive/VincentMoi/1ProjEV/StatModeling.v.01/fct_supr.R', chdir = TRUE)
source('~/Google Drive/VincentMoi/1ProjEV/StatModeling.v.01/fct_Unva.R', chdir = TRUE)
source('~/Google Drive/VincentMoi/1ProjEV/StatModeling.v.01/ORRR.R', chdir = TRUE)
source('~/Google Drive/VincentMoi/1ProjEV/StatModeling.v.01/Rank_test.R', chdir = TRUE)
source('~/Google Drive/VincentMoi/1ProjEV/StatModeling.v.01/Resum.cox.R', chdir = TRUE)
source('~/Google Drive/VincentMoi/1ProjEV/StatModeling.v.01/Resum.mcox.R', chdir = TRUE)
source('~/Google Drive/VincentMoi/1ProjEV/StatModeling.v.01/tab_cont.R', chdir = TRUE)


# noms des variables discretes
nomv<-c("BRC1","pTNM","pTNM1234","SeminalVesicleInvasion01","Margin01","StatusCRPC01","BCRType",
"RadiationType","Gleason_Sum","limphNodinv","CapsulP","meta","DeathOverall01")
newnam<-c("Biochemical Relapse","pTNM","pTNM1234","Seminal Gland Invasion","Margin status","StatusCRPC01",
"Biochemical Relapse Type","Radiation type","Gleason Sum","Lymph Node Invasion","Capsular Penetration",
"Bone Metastasis","Death")

Rst1<-Desc_fct(nomv,datmarq,namv1=newnam)
tabd1<-round(Rst1$tabR,digits=1)
require(Hmisc)
W<-latex(tabd1,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex')

# noms des variables continues
varc<-c("MeanBKi67","MeanTki67","MeanBp271","MeanTp271","MeanBp272","MeanTp272","datBRC",
"PSA","Age","datMeta","datsurv")
newvarc<-c("Ki67 normal ","ki67 tumoral","P27 normal (1)","P27 tumoral (1)","P27 normal (2)","P27 tumoral (2)","BRC date",
"PSA","Age","Metastasis date","Follow-up")

ll1<-Desc_Qte(varc,datmarq,namv=newvarc)

dtr1<-cbind(ll1$tabM,ll1$TabMed)
dtr2<-round(dtr1,digits=2)
W<-latex(dtr2,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)


# le programme qui permet de trouver le cutoff oppti

opicutoff<-function(mar1,stat){
    # marq : la valeur du marqueur
    # statu : la 1 ou 0
    # la fonction de l'optimisaation
    
    specif<-function(x){
        speci<-NULL
        for(xx in x){
            tab1<-table(rr,stat)
            rr<-ifelse(mar1>xx,1,0)+ifelse(is.na(mar1)==TRUE,NA,0)
            if(sum(rr)/length(rr[!is.na(rr)])==1){speci<-c(speci,table(stat)[2]/sum(tab1))}
            if(sum(rr)/length(rr[!is.na(rr)])==0){speci<-c(speci,table(stat)[1]/sum(tab1))}
            if(sum(rr)/length(rr[!is.na(rr)])>0 & sum(rr)/length(rr[!is.na(rr)])<1){
                speci<-c(speci,tab1[1,1]/sum(tab1)+tab1[2,2]/sum(tab1))}
        }
        return(speci)
    }
  
  return(specif)
}

# specificité et sensibilité
ff_SeSp<-function(mar1,stat){
    # marq : la valeur du marqueur
    # statu : la 1 ou 0
    # la fonction de l'optimisaation
    
    sen<-function(x){
        speci<-NULL
        for(xx in x){
            rr<-ifelse(mar1>=xx,1,0)+ifelse(is.na(mar1)==TRUE,NA,0)
             tab1<-table(rr,stat)
            if(sum(rr)/length(rr[!is.na(rr)])==1){speci<-c(speci,1)}
            if(sum(rr)/length(rr[!is.na(rr)])==0){speci<-c(speci,0)}
            if(sum(rr)/length(rr[!is.na(rr)])>0 & sum(rr)/length(rr[!is.na(rr)])<1){
                speci<-c(speci,tab1[2,2]/sum(tab1[,2]))}
        }
        return(speci)
    }
    
    specif<-function(x){
        speci<-NULL
        for(xx in x){
            rr<-ifelse(mar1>=xx,1,0)+ifelse(is.na(mar1)==TRUE,NA,0)
              tab1<-table(rr,stat)
            if(sum(rr)/length(rr[!is.na(rr)])==1){speci<-c(speci,0)}
            if(sum(rr)/length(rr[!is.na(rr)])==0){speci<-c(speci,1)}
            if(sum(rr)/length(rr[!is.na(rr)])>0 & sum(rr)/length(rr[!is.na(rr)])<1){
                speci<-c(speci,tab1[1,1]/sum(tab1[,1]))}
        }
        return(speci)
    }
    
    return(list(sensib=sen,specif=specif))
}

# Verification des programme
rr<-runif(250)
rrp<-rr
sta<-ifelse(rr>0.315158817,1,0)
ffe<-opicutoff(rrp,sta)
a=max(rrp)
curve(ffe(x),xlim=c(0,a))
abline(v=0.315158817,col=4)
abline(h=1,col=2)

oo<-ff_SeSp(rrp,sta)
se<-oo$sensib
sp<-oo$specif
curve(se(x),xlim=c(0,1),ylim=c(0,1))
curve(sp(x),add=TRUE,xlim=c(0,1),ylim=c(0,1),col=2)

cu<-seq(0,1,0.1)
sex<-se(cu)
spx<-1-sp(cu)
plot(sex~spx,typ='l')

# les fonction de transformation des variable

finv<-Vectorize(function(x){if(x<1){rr<-1}else{rr<-1/x}
return(rr)},'x')
flog<-Vectorize(function(x){if(x<1){rr<-0}else{rr<-log(x)}
    return(rr)},'x')
fsqrt<-function(x){sqrt(x)}
fid<-function(x){x}


# Recherche des meilleur cutoff
library(survival)
fcttm<-function(marq,stat,timv,ft,datmarq,vect,newvect,naout,vacox){
    naout1=paste(naout,"_Bi",sep="")
    naout2=paste(naout,"_Co",sep="")
    marmd<-ft(datmarq[[marq]])
    stat1<-datmarq[[stat]]
    a=max(marmd)
    rrt<-ff_SeSp(marmd,stat1)
    Sek<-rrt$sensib
    Spk<-rrt$specif
    ffcu<-function(){
        curve(Sek(x),xlim=c(0,a),ylim=c(0,1))
        curve(Spk(x),add=TRUE,xlim=c(0,a),ylim=c(0,1),col=2)
    }
    fuctM<-function(cutf,titr){
        datmarq[["maqd"]]<-factor(ifelse(marmd>=cutf,1,0),labels=c("Negative","Positive"))
        datmarq[[naout1]]<-factor(ifelse(marmd>=cutf,1,0))
        datmarq[[naout2]]<-marmd
        vacox1b=c(naout1,vacox)
        vacox1c=c(naout2,vacox)
        # la table descriptive
        ll<-Desc_fct(vect,datmarq,outcom1=naout1,namv1=newvect)
        #modèle de cox en discret
        tabb<-Fct_Rscox(vacox1b,vacox1b,stat,timv,datmarq)
        #modèle de cox en continue
        tabc<-Fct_Rscox(vacox1c,vacox1c,stat,timv,datmarq)
        # codage 1
        #Rank_test<-function(tim,even,ssta,data1,coln=NULL,titre1=NULL)
        llo<-Rank_test(timv,stat,"maqd",datmarq,titre1=titr)
        ffl<-llo$fup
        tdsc=ll$tabR
        return(list(flr=ffl,tabk=tdsc,tabb=tabb,tabc=tabc))
    }
    return(list(marmd=marmd,Sek=Sek,Spk=Spk,ffcu=ffcu,fuctM=fuctM))
    
}

# restriction des données
  datmarq$datBRC<-ifelse(datmarq$datBRC==0,3,datmarq$datBRC)
# datmarq<-datmarq[datmarq$BCRType!="FailedRP",]
  datmarq2<-datmarq

  vect=c("Gleason_Sum","pTNM1234","limphNodinv","CapsulP","SeminalVesicleInvasion01","Margin01",
"StatusCRPC01","BRC1","meta")
  newvect<-c("Gleason Sum","pTNM1234","Lymph Node Invasion","Capsular Penetration",
"Seminal Gland Invasion","Margin status","Status CRPC","Biochemical Relapse","Bone Metastasis")
  vacox=c("Age","pTNM1234","Gleason_Sum","Margin01")

# les arguments du modèle
  stat="BRC";timv="datBRC";naout="KI67.N";
#marq="MeanBKi67";;datmarq=datmarq2
#ft=fsqrt

# traitement du marqueur ki67 pas dans la tumeur (tissue normal)
llki67<-fcttm("MeanBKi67",stat,timv,fsqrt,datmarq2,vect,newvect,naout,vacox)
maqt<-llki67$marmd
fct1<-llki67$ffcu
fct1()

se1<-llki67$Sek
sp1<-llki67$Spk

fct2<-llki67$fuctM
ffrg<-fct2(1.007,"Ki67_N")
fctg<-ffrg$flr
tabd<-ffrg$tabk
w1<-latex(tabd,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)
tabb<-ffrg$tabb
w1<-latex(tabb,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)
tabc<-ffrg$tabc
w1<-latex(tabc,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)

ch1<-c("~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/grki67N.png")
png(file =ch1, bg = "transparent")
fctg(155,0.2,c("Time (Months)","BRC survival"))
dev.off()


# le ki67 dans la tumeur "finv,flog,fsqrt,fid"
library(xlsx)
ch1<-c("/Users/molierenguilemakao/tabRCox.xlsx")


stat="BRC";timv="datBRC";naout="KI67.T";
#marq="MeanBKi67";ft=fsqrt;datmarq=datmarq2

Rst1<-Desc_fct(c("Margin01","pTNM1234"),datmarq,"Gleason_Sum")
Rst2<-Desc_fct(c("Gleason_Sum","pTNM1234"),datmarq,"Margin01")

# traitement du marqueur ki67 pas dans la tumeur
vacox=c("Age","Margin01","pTNM1234","PSA")# "Gleason_Sum",
vacox=c("Age","Gleason_Sum","PSA")# "Gleason_Sum",
llki67<-fcttm("MeanTki67",stat,timv,fsqrt,datmarq2,vect,newvect,naout,vacox)
maqt<-llki67$marmd
fct1<-llki67$ffcu
fct1()

se1<-llki67$Sek
sp1<-llki67$Spk

fct2<-llki67$fuctM
#ffrg<-fct2(1.36,"Ki67_T")
ffrg<-fct2(2.3,"Ki67_T")
fctg<-ffrg$flr
tabd<-ffrg$tabk
w1<-latex(tabd,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)
tabb<-ffrg$tabb
write.xlsx(tabb,file=ch1,sheetName="Ancoxk67d")
w1<-latex(tabb,file='~/Google Drive/VincentMoi/1ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/exmp.tex')
tabc<-ffrg$tabc
write.xlsx(tabc,file=ch1,sheetName="Ancoxk67c",append=TRUE)

w1<-latex(tabc,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)

ch1<-c("~/Google Drive/VincentMoi/1ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/grki67T.png")
png(file =ch1, bg = "transparent")
fctg(155,0.2,c("Time (Months)","BRC survival"))
dev.off()

# le marqueur p27 hors de la tumeur "finv,flog,fsqrt,fid"

stat="BRC";timv="datBRC";naout="p27.N";
#marq="MeanBKi67";ft=fsqrt;datmarq=datmarq2
# traitement du marqueur ki67 pas dans la tumeur

vacox=c("Age","pTNM1234","Gleason_Sum","Margin01","limphNodinv","CapsulP","SeminalVesicleInvasion01")
llki67<-fcttm("MeanBp271",stat,timv,finv,datmarq2,vect,newvect,naout,vacox)
maqt<-llki67$marmd
fct1<-llki67$ffcu
fct1()

se1<-llki67$Sek
sp1<-llki67$Spk

fct2<-llki67$fuctM
ffrg<-fct2(0.0343,"P27_N")
fctg<-ffrg$flr
tabd<-ffrg$tabk
w1<-latex(tabd,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)
tabb<-ffrg$tabb
w1<-latex(tabb,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)
tabc<-ffrg$tabc
w1<-latex(tabc,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)

ch1<-c("~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/p27N.png")
png(file =ch1, bg = "transparent")
fctg(155,0.2,c("Time (Months)","BRC survival"))
dev.off()




stat="BRC";timv="datBRC";naout="p27.T";
#marq="MeanBKi67";ft=fsqrt;datmarq=datmarq2
# traitement du marqueur ki67 pas dans la tumeur
vacox=c("Age","Gleason_Sum","Margin01")
vacox=c("Age","Gleason_Sum","PSA")

llki67<-fcttm("MeanTp271",stat,timv,finv,datmarq2,vect,newvect,naout,vacox)
maqt<-llki67$marmd
fct1<-llki67$ffcu
fct1()

se1<-llki67$Sek
sp1<-llki67$Spk

fct2<-llki67$fuctM
ffrg<-fct2(0.0289,"P27_T")
fctg<-ffrg$flr
fctg(155,0.2,c("Time (Months)","BRC survival"))

tabd<-ffrg$tabk
w1<-latex(tabd,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)
tabb<-ffrg$tabb
write.xlsx(tabb,file=ch1,sheetName="Ancoxp27b",append=TRUE)
w1<-latex(tabb,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)
tabc<-ffrg$tabc
write.xlsx(tabc,file=ch1,sheetName="Ancoxp27c1",append=TRUE)
w1<-latex(tabc,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)

ch1<-c("~/Google Drive/VincentMoi/1ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/p27T.png")
png(file=ch1, bg = "transparent")
fctg(100,0.2,c("Time (Months)","BRC survival"))
dev.off()




# Pour le marqueur p27 hors de la tumeur (code 2)
vacox=c("Age","pTNM1234","Gleason_Sum","Margin01","limphNodinv","CapsulP","SeminalVesicleInvasion01")
stat="BRC";timv="datBRC";naout="p27.N";
#marq="MeanBKi67";ft=fsqrt;datmarq=datmarq2
# traitement du marqueur ki67 pas dans la tumeur
#vacox=c("Age","Gleason_Sum","Margin01")

llki67<-fcttm("MeanBp272",stat,timv,finv,datmarq2,vect,newvect,naout,vacox)
maqt<-llki67$marmd
fct1<-llki67$ffcu
fct1()

se1<-llki67$Sek
sp1<-llki67$Spk

fct2<-llki67$fuctM
ffrg<-fct2(0.02599,"P27_N")
fctg<-ffrg$flr
fctg(155,0.2,c("Time (Months)","BRC survival"))

tabd<-ffrg$tabk
w1<-latex(tabd,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)
tabb<-ffrg$tabb
w1<-latex(tabb,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)
tabc<-ffrg$tabc
w1<-latex(tabc,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)

ch1<-c("~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/p272N.png")
png(file =ch1, bg = "transparent")
fctg(100,0.2,c("Time (Months)","BRC survival"))
dev.off()







# la tumeur
# Pour le marqueur p27 hors de la tumeur (code 2)
#vacox=c("Age","pTNM1234","Gleason_Sum","Margin01","limphNodinv","CapsulP","SeminalVesicleInvasion01")
stat="BRC";timv="datBRC";naout="p27.T";
#marq="MeanBKi67";ft=fsqrt;datmarq=datmarq2
# traitement du marqueur ki67 pas dans la tumeur avec fsqrt 8.1
vacox=c("Age","Gleason_Sum","PSA")

llki67<-fcttm("MeanTp272",stat,timv,fsqrt,datmarq2,vect,newvect,naout,vacox)
maqt<-llki67$marmd
fct1<-llki67$ffcu
fct1()

se1<-llki67$Sek
sp1<-llki67$Spk

fct2<-llki67$fuctM
ffrg<-fct2(8.1,"P27_T")
#ffrg<-fct2(8.7,"P27_T")
fctg<-ffrg$flr
fctg(155,0.2,c("Time (Months)","BRC survival"))

tabd<-ffrg$tabk
w1<-latex(tabd,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)
tabb<-ffrg$tabb
w1<-latex(tabb,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)
tabc<-ffrg$tabc
w1<-latex(tabc,file='~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/exmp.tex',append=TRUE)

ch1<-c("~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/p272T1.png")
png(file =ch1, bg = "transparent")
fctg(100,0.2,c("Time (Months)","BRC survival"))
dev.off()

llki67<-fcttm("MeanBp272","BRC",finv,datmarq)
maqt<-llki67$marmd
fct1<-llki67$ffcu
fct1()

se1<-llki67$Sek
sp1<-llki67$Spk

fct2<-llki67$fuctM
ch1<-c("~/Google Drive/VincentMoi/ProjEV/Proj_CPCBN_Terry_fox/Dos_Redac/docR_tex/grki67t.png")
png(file =ch1, bg = "transparent")
fct2(0.02,15,0.4,"p27.T_cod1") #76-80 pour finv


# le test des donnéés
stat="BRC";timv="datBRC"
vec<-c("pTNM1234","Margin01","PSA","Age","Gleason_Sum")
rst<-Fct_Rscox(c("PSA","Age"),vec,stat,timv,datmarq)
write.xlsx(rst,file="/Users/molierenguilemakao/modcox.xlsx")













































