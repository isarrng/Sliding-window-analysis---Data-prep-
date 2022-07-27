
library(ggmap)
library(pracma)
library(snpStats)

anno=read.csv("v50.0_1240K_public1.anno",header=FALSE,sep="\t")
PLINKfam=read.table("v50.0_1240k_public_plink.fam")

fromBP=7000
toBP=0

#Input values for maximum and minimum latitude and longitude#

maxlat=60
minlat=50
maxlon=1.5
minlon=-10

#Input latitudes and longtitudes for exclusion#

exmaxlat=51
exminlat=50
exmaxlon=1.5
exminlon=-2.5

ex1maxlat=43
ex1minlat=38
ex1maxlon=10
ex1minlon=8

ex2maxlat=40
ex2minlat=38
ex2maxlon=4.5
ex2minlon=1.1


#Subsetting individuals. Addition logic statements added below if required#

anno$V16=as.numeric(anno$V16)
anno$V17=as.numeric(anno$V17)
individuals=subset(anno,anno$V16>minlat&anno$V16<maxlat&anno$V17>minlon&
                     anno$V17<maxlon&anno$V9>toBP&anno$V9<fromBP)

#Disable following 3 lines if nothing is to be excluded#

individuals=subset(individuals, !(individuals$V16>exminlat&individuals$V16<exmaxlat&individuals$V17>exminlon&
                                    individuals$V17<exmaxlon))

individuals=subset(individuals, !(individuals$V16>ex1minlat&individuals$V16<ex1maxlat&individuals$V17>ex1minlon&
                                    individuals$V17<ex1maxlon))

individuals=subset(individuals, !(individuals$V16>ex2minlat&individuals$V16<ex2maxlat&individuals$V17>ex2minlon&
                                    individuals$V17<ex2maxlon))


individualsGeoDate=data.frame(individuals$V2,individuals$V9,individuals$V10,
                              individuals$V16,individuals$V17)
names(individualsGeoDate)[names(individualsGeoDate)=="individuals.V2"]="IID"
names(individualsGeoDate)[names(individualsGeoDate)=="individuals.V9"]="meanBP"
names(individualsGeoDate)[names(individualsGeoDate)=="individuals.V10"]="sdBP"
names(individualsGeoDate)[names(individualsGeoDate)=="individuals.V16"]="lat"
names(individualsGeoDate)[names(individualsGeoDate)=="individuals.V17"]="lon"
individualsGeoDate$FID=PLINKfam$V1[match(individualsGeoDate$IID,PLINKfam$V2)]
individualsGeoDate=individualsGeoDate[, c(6,1,2,3,4,5)]

write.table(individualsGeoDate,file="individualsGeoDate.txt",col.names=TRUE,
            row.names=FALSE)


#Create map#

euromid=array(data=NA,dim=c(1,2,1))
euromid[1,1,1]=11.556853
euromid[1,2,1]=50.045074
euromid=c(lat=50.045074,lon=11.556853)
euromap=get_map(location=euromid,zoom=4,crop=FALSE)

ggmap(euromap)+
  geom_point(aes(x=lon,y=lat),data=individualsGeoDate,alpha=0.5,color="darkred",
             size=1)


#Creating file for PLINK#

individualsP=data.frame(individuals$V2)
names(individualsP)[names(individualsP)=="individuals.V2"]="V2"
individualsP$V1=PLINKfam$V1[match(individuals$V2,PLINKfam$V2)]
individualsP$V3=PLINKfam$V3[match(individuals$V2,PLINKfam$V2)]
individualsP$V4=PLINKfam$V4[match(individuals$V2,PLINKfam$V2)]
individualsP$V5=PLINKfam$V5[match(individuals$V2,PLINKfam$V2)]
individualsP$V6=PLINKfam$V6[match(individuals$V2,PLINKfam$V2)]
individualsP=individualsP[, c(2,1,3,4,5,6)]


write.table(individualsP,file="prePLINK_ind_7to0_UK.txt",col.names=TRUE,
            row.names=FALSE,quote=FALSE)


################################AFTER-PLINK####################################


#Input PLINK .fam .bed .bim files#

individualsGeoDate=read.table("individualsGeoDate.txt")
fam=read.table("ind_7to0_m45_g50.fam")
indGeno=read.plink("ind_7to0_m45_g50.bed",
                   "ind_7to0_m45_g50.bim",
                   "ind_7to0_m45_g50.fam")

#Create temporal distance matrix#

fam$meanBP=individualsGeoDate$meanBP[match(fam$V2,individualsGeoDate$IID)]
fam$sdBP=individualsGeoDate$sdBP[match(fam$V2,individualsGeoDate$IID)]
fam$fromBP=fam$meanBP+fam$sdBP
fam$toBP=fam$meanBP-fam$sdBP


indDates=data.frame(fam$V2,fam$fromBP,fam$toBP)
names(indDates)[names(indDates)=="fam.V2"]="IID"
names(indDates)[names(indDates)=="fam.fromBP"]="fromBP"
names(indDates)[names(indDates)=="fam.toBP"]="toBP"

write.table(indDates,file="ind_7to0_UK_m45_g50.txt",col.names=TRUE,
            row.names=FALSE)


#Create geographic distance matrix#

fam$lat=individualsGeoDate$lat[match(fam$V2,individualsGeoDate$IID)]
fam$lon=individualsGeoDate$lon[match(fam$V2,individualsGeoDate$IID)]


ukmid=c(lat=54,lon=-4.7)
ukmap=get_map(location=ukmid,zoom=5,crop=FALSE)

ggmap(ukmap)+
  geom_point(aes(x=lon,y=lat),data=fam,alpha=0.5,color="darkred",
             size=1)




indGeo=data.frame(fam$V2,fam$lat,fam$lon)

names(indGeo)[names(indGeo)=="fam.V2"]="IID"
names(indGeo)[names(indGeo)=="fam.lat"]="lat"
names(indGeo)[names(indGeo)=="fam.lon"]="lon"

indGeoMatrix=as.data.frame(matrix(ncol=nrow(indGeo),nrow=nrow(indGeo)))

colnames(indGeoMatrix)=c(indGeo$IID)


for(i in 1:nrow(indGeo)) {
  for(j in 1:nrow(indGeo)) {
    indGeoMatrix[i,j]=
      haversine(loc1=as.numeric(c(indGeo[i,2],indGeo[i,3])),
                loc2=as.numeric(c(indGeo[j,2],indGeo[j,3])))
  }
}

write.table(indGeoMatrix,file="ind_7to0_UK_m45_g50_Geo.txt",col.names=TRUE,
            row.names=FALSE)


#Create genetic distance matrix#

indGeno_ibs=ibsCount(indGeno$genotypes,
                     uncertain=FALSE)

indGeno_ibs_dist=ibsDist(indGeno_ibs)

indGeno_ibs_dist_m=as.matrix(indGeno_ibs_dist)

indGeno_ibs_dist_df=as.data.frame(indGeno_ibs_dist_m)

write.table(indGeno_ibs_dist_df,file="ind_7to0_UK_m45_g50_Gen.txt",col.names=TRUE,
            row.names=FALSE)

