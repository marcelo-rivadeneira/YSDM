#' LazySDM
#'
#' Black-box correlational SDM, with minimum specification steps. Calibrates SDM in time1 and projects onto time2.
#' @param dato data.frame containing geographic coordinates of a species presence
#' @param buff geographic buffer (in km) where to generate pseudo-absence. It is also the area where the model is tested.
#' @param stck1 Raster stack (terra class) with environmental variables for time 1
#' @param stck2 Raster stack (terra class) with environmental variables for time 2
#'
#' @return Four objects
#'      summary: Main diagnostic statistics
#'      target.area: raster of the buffer zone
#'      occupancy.time1: raster of estimated habitability for time 1
#'      occupancy.time2: raster of estimated habitability for time 2
#' @import terra
#' @import sp
#' @import blockCV
#' @import SDMtune
#' @import ENMeval
#' @import usdm
#' @import PresenceAbsence
#' @import modEvA
#' @export
#'
#' @examples
#' lazySDM(whiteshark,buff=100,stck1,stck2)
#'
lazySDM=function(dato,buff=500,stck1,stck2)
{
  options(warn=-1)
  x=dato$x
  y=dato$y
  coords=data.frame(x,y)
  coordinates(coords)= ~x + y
  crs(coords)="+proj=longlat +datum=WGS84"
  buf=vect(buffer(coords,(buff*1e3))) # default, 500 km buffer
  mask=stck1[[1]]  # Create a mask with target resolution and extent from climate layers
  values(mask)[!is.na(values(mask))]=1
  maskedbuf=mask(mask,buf) # Set all raster cells outside the buffer to NA.
  ext.coods=intersect(ext(coords), ext(maskedbuf))
  maskedbuf2=crop(maskedbuf, ext.coods)
  bg_dat=spatSample(maskedbuf2,size=length(x),na.rm=T,as.points=T)
  bg_dat=crds(bg_dat)
  pr_dat=data.frame(x,y)
  puntos=rbind(pr_dat,bg_dat)

  # Autocorrelation range of environmental variables
  sac=cv_spatial_autocor(crop(stck1,maskedbuf),num_sample=1000,progress=T,plot=F)
  msac=median(sac$range_table$range)/1000

  ## Prepare the tunning
  totuneo=prepareSWD(species="bla",p = pr_dat,a = bg_dat, env = stck1)

  ## spatial folding (4 checkerboard) for cross-validation
  check_folds=get.checkerboard2(occs = pr_dat,env=stck1,bg = bg_dat,aggregation.factor = 4)
  folds.p=check_folds$occs.grp
  block_folds_formatted1=matrix(ncol=length(unique(folds.p)),nrow=length(folds.p),F)
  block_folds_formatted2=matrix(ncol=length(unique(folds.p)),nrow=length(folds.p),T)
  for(i in 1:nrow(block_folds_formatted1))
  {
    block_folds_formatted1[i,folds.p[i]]=T
    block_folds_formatted2[i,folds.p[i]]=F
  }
  block_folds_formatted=list(block_folds_formatted1,block_folds_formatted2)
  names(block_folds_formatted)=c("train","test")

  # Cross-validated trained RF model with hyperparameter tunning
  randfo.model=train("RF", data = totuneo, folds = block_folds_formatted,verbose=F)
  randfo.auc.test=SDMtune::auc(randfo.model,test=T)


  # Thresholding
  rf.predicted=predict(randfo.model, data = totuneo)
  spatcoord=totuneo@coords
  coordinates(spatcoord)=~X+Y
  umbral=cbind(1:length(totuneo@pa),sobs=totuneo@pa,
               avg.prob=rf.predicted)
  av.thr=optimal.thresholds(umbral)[3,2]

  # variable importance
  impset=data.frame(varname=names(stck1))
  impset$numname=1:nrow(impset)

  vimp=SDMtune::varImp(randfo.model,permut=10)
  #vimp$numname=impset$numname[match(vimp$Variable,impset$varname)]
  #vimp=vimp[order(vimp$numname),]

  # temporal collinearity shift between predictors
  messy.stck1=data.frame(totuneo@data)
  messy.stck2=terra::extract(x=stck2,y=vect(spatcoord,crs="+proj=longlat +datum=WGS84"),ID=F)
  vif.stck1=stats::cor(messy.stck1)
  vif.stck2=cor(messy.stck2)
  cor.present.stck1=cor(as.dist(vif.stck1),as.dist(vif.stck2))

  # VIF of predictors
  vif1=vif(messy.stck1)
  vif2=vif(messy.stck2)

  v.vif1=vif1$Variables[which(vif1$VIF==max(vif1$VIF))]
  l.vif1=vif1$VIF[which(vif1$VIF==max(vif1$VIF))]

  v.vif2=vif2$Variables[which(vif2$VIF==max(vif2$VIF))]
  l.vif2=vif2$VIF[which(vif2$VIF==max(vif2$VIF))]


  #MESS
  mess.p1=modEvA::MESS(V = totuneo@data[totuneo@pa==1,], P = totuneo@data[totuneo@pa==0,],verbosity=0)
  mess1=sum(ifelse(mess.p1$TOTAL<0,1,0))/nrow(mess.p1)

  mess.p2=modEvA::MESS(V = totuneo@data[totuneo@pa==1,], P = messy.stck2,verbosity=0)
  mess2=sum(ifelse(mess.p2$TOTAL<0,1,0))/nrow(mess.p2)


  ## Predicted distribution in the target area for present and future scenarios
  hab.stck1=predict(randfo.model, data = stck1)
  hab.stck1[hab.stck1>av.thr]=1;hab.stck1[hab.stck1<1]=0

  hab.stck2=predict(randfo.model, data = stck2)
  hab.stck2[hab.stck2>av.thr]=1;hab.stck2[hab.stck2<1]=0


  ## Predicted area of occupancy (AOO) in the Southern Ocean for each scenario
  aoo.stck1=global(cellSize(hab.stck1,unit="km")*hab.stck1, "sum",na.rm=T)$sum
  aoo.stck2=global(cellSize(hab.stck2,unit="km")*hab.stck2, "sum",na.rm=T)$sum


  outta=data.frame(nobs=nrow(pr_dat), buff,auc.test=randfo.auc.test,
                   topvarimp=vimp$Variable[1],topimp=vimp$Permutation_importance[1],msac,
                   v.vif1,l.vif1,v.vif2,l.vif1, cor.present.stck1,
                   mess1, mess2, aoo.t1=aoo.stck1,aoo.t2=aoo.stck2)
  outta=data.frame(t(outta))

  salida=list(outta,maskedbuf,hab.stck1,hab.stck2)
  names(salida)=c("summary","target.area","occupancy.time1","occupancy.time2")
  return(salida)
  options(warn=0)
}
