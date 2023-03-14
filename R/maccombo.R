# Maxcombo Functions

#'
#' @importFrom magrittr '%>%'
#'
oogetdataframemsdataprocessedvRECODED=function(oodataframe,oostringorsymbolid="id",oostringorsymboltreated="treated",oostringorsymbolAtime="Atime",oostringorsymbolBtime="Btime",oostringorsymbolBobserved="Bobserved",oostringorsymbolCtime="Ctime",oostringorsymbolCobserved="Cobserved")
{
  ooexpressionid=rlang::ensym(oostringorsymbolid)
  ooexpressiontreated=rlang::ensym(oostringorsymboltreated)
  ooexpressionAtime=rlang::ensym(oostringorsymbolAtime)
  ooexpressionBtime=rlang::ensym(oostringorsymbolBtime)
  ooexpressionBobserved=rlang::ensym(oostringorsymbolBobserved)
  ooexpressionCtime=rlang::ensym(oostringorsymbolCtime)
  ooexpressionCobserved=rlang::ensym(oostringorsymbolCobserved)

  oomatrixinttmat=base::matrix(NA_integer_,nrow=3L,ncol=3L)
  oomatrixinttmat[[1L,2L]]=1L; oomatrixinttmat[[1L,3L]]=2L;
  base::dimnames(oomatrixinttmat)=base::list(from=c("A","B","C"),to=c("A","B","C") )

  oodataframesupport=oodataframe %>% dplyr::select(-c(!!ooexpressionAtime,!!ooexpressionBtime,!!ooexpressionBobserved,!!ooexpressionCtime,!!ooexpressionCobserved) ) #auxiliary covariates, since mstate::msprep keep= argument seems to be broken
  oodataframemsdata=mstate::msprep(time=c(rlang::as_string(ooexpressionAtime),rlang::as_string(ooexpressionBtime),rlang::as_string(ooexpressionCtime) ),
                                   status=c(NA,rlang::as_string(ooexpressionBobserved),rlang::as_string(ooexpressionCobserved) ),
                                   data=oodataframe,
                                   trans=oomatrixinttmat,
                                   id=rlang::as_string(ooexpressionid),
                                   start=base::list(state=base::rep(1L,base::nrow(oodataframe) ),time=oodataframe[[rlang::as_string(ooexpressionAtime)]] )
  )
  oodataframemsdata=base::merge(oodataframemsdata,oodataframesupport,by=rlang::as_string(ooexpressionid) )
  oodataframemsdata
}


oogetresultsbasicv2vRECODED=function(oodataframe,oofunctionweightasafunctionofstminus=function(oodoublestminus){ base::return(1.0) },oostringorsymbolid="id",oostringorsymboltreated="treated",oostringorsymbolAtime="Atime",oostringorsymbolBtime="Btime",oostringorsymbolBobserved="Bobserved",oostringorsymbolCtime="Ctime",oostringorsymbolCobserved="Cobserved")
{
  ooexpressionid=rlang::ensym(oostringorsymbolid)
  ooexpressiontreated=rlang::ensym(oostringorsymboltreated)
  ooexpressionAtime=rlang::ensym(oostringorsymbolAtime)
  ooexpressionBtime=rlang::ensym(oostringorsymbolBtime)
  ooexpressionBobserved=rlang::ensym(oostringorsymbolBobserved)
  ooexpressionCtime=rlang::ensym(oostringorsymbolCtime)
  ooexpressionCobserved=rlang::ensym(oostringorsymbolCobserved)

  oodataframetrans1=oodataframe %>% oogetdataframemsdataprocessedvRECODED(oostringorsymbolid = !!ooexpressionid,oostringorsymboltreated = !!ooexpressiontreated,oostringorsymbolAtime = !!ooexpressionAtime,oostringorsymbolBtime = !!ooexpressionBtime,oostringorsymbolBobserved = !!ooexpressionBobserved,oostringorsymbolCtime = !!ooexpressionCtime,oostringorsymbolCobserved = !!ooexpressionCobserved) %>% dplyr::filter(.data$trans==1L) #see explanation of why I am using .data$trans==1L instead of trans==1L here in the following guide about how to use dplyr from within a package while avoiding R CMD CHECK notes: https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html
  oovecinttreated=oodataframetrans1[[rlang::as_string(ooexpressiontreated) ]]
  oovecdoublemintc=oodataframetrans1$time
  oovecboolobservedt=oodataframetrans1$status

  ooSurv=survival::Surv(time=oovecdoublemintc,event=oovecboolobservedt,type="right")
  oosurvfit=survival::survfit(formula=ooSurv~1)
  #base::unclass(oosurvfit)

  #this particular variable is used in many of the weight specifications, so we calculate it here
  oovecdoublesurvivaltminus=c(1.0,oosurvfit$surv)[1L:(base::length(oosurvfit$surv) ) ]

  #oovecinty1plusy2=oosurvfit$n.risk
  #oosurvfit$time %>% purrr::map_int(~base::sum(oovecdoublemintc >= .x)  ) #should equal the above
  oovecinty1=oosurvfit$time %>% purrr::map_int(~base::sum(oovecdoublemintc[oovecinttreated==0L] >= .x)  )
  oovecinty2=oosurvfit$time %>% purrr::map_int(~base::sum(oovecdoublemintc[oovecinttreated==1L] >= .x)  )

  #oovecintdn1plusdn2=oosurvfit$n.event
  #oosurvfit$time %>% purrr::map_int(~base::sum(oovecdoublemintc==.x & oovecboolobservedt==TRUE) ) #should equal the above
  oovecintdn1=oosurvfit$time %>% purrr::map_int(~base::sum(oovecdoublemintc[oovecinttreated==0L]==.x & oovecboolobservedt[oovecinttreated==0L]==TRUE) )
  oovecintdn2=oosurvfit$time %>% purrr::map_int(~base::sum(oovecdoublemintc[oovecinttreated==1L]==.x & oovecboolobservedt[oovecinttreated==1L]==TRUE) )

  oodataframebytimepoint=dplyr::tibble(oovecdoubletime=oosurvfit$time,oovecdoublesurvivaltminus,oovecinty1,oovecinty2,oovecintdn1,oovecintdn2)
  oodataframebytimepointreduced=oodataframebytimepoint %>% dplyr::filter(oovecinty1 > 0L & oovecinty2 > 0L & (oovecintdn1 > 0L | oovecintdn2 > 0L) )

  oovecdoubleweights=oodataframebytimepointreduced$oovecdoublesurvivaltminus %>% purrr::map_dbl(oofunctionweightasafunctionofstminus)

  #see "Lee 2007 On the versatility of the combination of weighted log-rank statistics.pdf"
  oodoublestatisticnotnormalized=base::sum( oovecdoubleweights * (oodataframebytimepointreduced$oovecinty1*oodataframebytimepointreduced$oovecinty2)/(oodataframebytimepointreduced$oovecinty1 + oodataframebytimepointreduced$oovecinty2)*(oodataframebytimepointreduced$oovecintdn1/oodataframebytimepointreduced$oovecinty1 - oodataframebytimepointreduced$oovecintdn2/oodataframebytimepointreduced$oovecinty2) )

  #see "Lee 2007 On the versatility of the combination of weighted log-rank statistics.pdf"
  oodoublevarianceestimateofunnormalizedstatistic=base::sum( oovecdoubleweights * oovecdoubleweights * (oodataframebytimepointreduced$oovecinty1*oodataframebytimepointreduced$oovecinty2)/(oodataframebytimepointreduced$oovecinty1 + oodataframebytimepointreduced$oovecinty2)*(1 - (oodataframebytimepointreduced$oovecintdn1 + oodataframebytimepointreduced$oovecintdn2 - 1)/(oodataframebytimepointreduced$oovecinty1 + oodataframebytimepointreduced$oovecinty2 - 1) ) * (oodataframebytimepointreduced$oovecintdn1 + oodataframebytimepointreduced$oovecintdn2)/(oodataframebytimepointreduced$oovecinty1 + oodataframebytimepointreduced$oovecinty2) )

  oolistreturnvalue=base::list(oodoublestatisticnotnormalized=oodoublestatisticnotnormalized,oodoublevarianceestimateofunnormalizedstatistic=oodoublevarianceestimateofunnormalizedstatistic,oodoublestatisticnormalized=oodoublestatisticnotnormalized/base::sqrt(oodoublevarianceestimateofunnormalizedstatistic) )

  oolistreturnvalue

}


oogetdoublemaxcomboteststatistic=function(oodataframe,oolistfunctionweightasafunctionofstminus=base::list(function(oodoublestminus){ base::return(1.0) } ),oostringorsymbolid="id",oostringorsymboltreated="treated",oostringorsymbolAtime="Atime",oostringorsymbolBtime="Btime",oostringorsymbolBobserved="Bobserved",oostringorsymbolCtime="Ctime",oostringorsymbolCobserved="Cobserved")

{

  ooexpressionid=rlang::ensym(oostringorsymbolid)

  ooexpressiontreated=rlang::ensym(oostringorsymboltreated)

  ooexpressionAtime=rlang::ensym(oostringorsymbolAtime)

  ooexpressionBtime=rlang::ensym(oostringorsymbolBtime)

  ooexpressionBobserved=rlang::ensym(oostringorsymbolBobserved)

  ooexpressionCtime=rlang::ensym(oostringorsymbolCtime)

  ooexpressionCobserved=rlang::ensym(oostringorsymbolCobserved)



  if(base::length(oolistfunctionweightasafunctionofstminus) < 1L)

  {

    base::stop("oolistfunctionweightasafunctionofstminus must have length greater than or equal to 1.")

  }



  oovecdoublestandardizedweightedlogrankstatistics=oolistfunctionweightasafunctionofstminus %>% purrr::map_dbl(~oogetresultsbasicv2vRECODED(oodataframe=oodataframe,oofunctionweightasafunctionofstminus=.x,oostringorsymbolid = !!ooexpressionid,oostringorsymboltreated = !!ooexpressiontreated,oostringorsymbolAtime = !!ooexpressionAtime,oostringorsymbolBtime = !!ooexpressionBtime,oostringorsymbolBobserved = !!ooexpressionBobserved,oostringorsymbolCtime = !!ooexpressionCtime,oostringorsymbolCobserved = !!ooexpressionCobserved)[["oodoublestatisticnormalized"]] )


  oovecdoublestandardizedweightedlogrankstatistics = abs(oovecdoublestandardizedweightedlogrankstatistics)

  base::return(oovecdoublestandardizedweightedlogrankstatistics %>% base::max( ) )

}



oogetresultsbasiccovariance3vRECODED=function(oodataframet1,oodataframet2,oofunctionweightasafunctionofstminus1=function(stminus){ base::return(1.0) },oofunctionweightasafunctionofstminus2=function(stminus){ base::return(1.0) }, oostringorsymbolid="id",oostringorsymboltreated="treated",oostringorsymbolAtime="Atime",oostringorsymbolBtime="Btime",oostringorsymbolBobserved="Bobserved",oostringorsymbolCtime="Ctime",oostringorsymbolCobserved="Cobserved")
{
  ooexpressionid=rlang::ensym(oostringorsymbolid)
  ooexpressiontreated=rlang::ensym(oostringorsymboltreated)
  ooexpressionAtime=rlang::ensym(oostringorsymbolAtime)
  ooexpressionBtime=rlang::ensym(oostringorsymbolBtime)
  ooexpressionBobserved=rlang::ensym(oostringorsymbolBobserved)
  ooexpressionCtime=rlang::ensym(oostringorsymbolCtime)
  ooexpressionCobserved=rlang::ensym(oostringorsymbolCobserved)


  oodataframetrans1=oodataframet1 %>% oogetdataframemsdataprocessedvRECODED(oostringorsymbolid = !!ooexpressionid,oostringorsymboltreated = !!ooexpressiontreated,oostringorsymbolAtime = !!ooexpressionAtime,oostringorsymbolBtime = !!ooexpressionBtime,oostringorsymbolBobserved = !!ooexpressionBobserved,oostringorsymbolCtime = !!ooexpressionCtime,oostringorsymbolCobserved = !!ooexpressionCobserved) %>% dplyr::filter(.data$trans==1L) #see explanation of why I am using .data$trans==1L instead of trans==1L here in the following guide about how to use dplyr from within a package while avoiding R CMD CHECK notes: https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html
  oovecinttreated=oodataframetrans1[[rlang::as_string(ooexpressiontreated) ]]
  oovecdoublemintc=oodataframetrans1$time
  oovecboolobservedt=oodataframetrans1$status

  ooSurv=survival::Surv(time=oovecdoublemintc,event=oovecboolobservedt,type="right")
  oosurvfit=survival::survfit(formula=ooSurv~1)
  #base::unclass(oosurvfit)

  #this particular variable is used in many of the weight specifications, so we calculate it here
  oovecdoublesurvivaltminus=c(1.0,oosurvfit$surv)[1:(base::length(oosurvfit$surv) ) ]

  #oovecinty1plusy2=oosurvfit$n.risk
  #oosurvfit$time %>% purrr::map_int(~base::sum(oovecdoublemintc >= .x)  ) #should equal the above
  oovecinty1=oosurvfit$time %>% purrr::map_int(~base::sum(oovecdoublemintc[oovecinttreated==0] >= .x)  )
  oovecinty2=oosurvfit$time %>% purrr::map_int(~base::sum(oovecdoublemintc[oovecinttreated==1] >= .x)  )

  #oovecintdn1plusdn2=oosurvfit$n.event
  #oosurvfit$time %>% purrr::map_int(~base::sum(oovecdoublemintc==.x & oovecboolobservedt==TRUE) ) #should equal the above
  oovecintdn1=oosurvfit$time %>% purrr::map_int(~base::sum(oovecdoublemintc[oovecinttreated==0]==.x & oovecboolobservedt[oovecinttreated==0]==TRUE) )
  oovecintdn2=oosurvfit$time %>% purrr::map_int(~base::sum(oovecdoublemintc[oovecinttreated==1]==.x & oovecboolobservedt[oovecinttreated==1]==TRUE) )

  oodataframebytimepoint=dplyr::tibble(oovecdoubletime=oosurvfit$time,oovecdoublesurvivaltminus,oovecinty1,oovecinty2,oovecintdn1,oovecintdn2)
  oodataframebytimepointreduced=oodataframebytimepoint %>% dplyr::filter(oovecinty1 > 0 & oovecinty2 > 0 & (oovecintdn1 > 0 | oovecintdn2 > 0) )

  oovecdoubleweights=oodataframebytimepointreduced$oovecdoublesurvivaltminus %>% purrr::map_dbl(oofunctionweightasafunctionofstminus1)


  #### This below part is new.  Now we also need the weight function from the second timepoint. ####

  oodataframetrans1t2=oodataframet2 %>% oogetdataframemsdataprocessedvRECODED(oostringorsymbolid = !!ooexpressionid,oostringorsymboltreated = !!ooexpressiontreated,oostringorsymbolAtime = !!ooexpressionAtime,oostringorsymbolBtime = !!ooexpressionBtime,oostringorsymbolBobserved = !!ooexpressionBobserved,oostringorsymbolCtime = !!ooexpressionCtime,oostringorsymbolCobserved = !!ooexpressionCobserved) %>% dplyr::filter(.data$trans==1L) #see explanation of why I am using .data$trans==1L instead of trans==1L here in the following guide about how to use dplyr from within a package while avoiding R CMD CHECK notes: https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html
  oosurvfitt2=survival::survfit(formula=survival::Surv(time=oodataframetrans1t2$time,event=oodataframetrans1t2$status,type="right")~1)

  oovecdoublesurvivaltminust2=c(1.0,oosurvfitt2$surv)[1:(base::length(oosurvfitt2$surv) ) ]

  oodataframeworkingt2=dplyr::tibble(time=oosurvfitt2$time,oovecdoublesurvivalminusATt2=oovecdoublesurvivaltminust2)

  oogetdoublefromdoublesurvivalminusAtt2=function(oodoubletime)
  {
    ooointcolumnindex=base::match(oodoubletime,oodataframeworkingt2$time)
    if(!base::is.na(ooointcolumnindex) )
    {
      base::return(oodataframeworkingt2$oovecdoublesurvivalminusATt2[[ooointcolumnindex]] )
    }else
    {
      oointcolumnindex=base::which.min( base::abs(oodataframeworkingt2$time - oodoubletime) )
      base::return(oodataframeworkingt2$oovecdoublesurvivalminusATt2[[oointcolumnindex]] )
    }
  }

  oovecdoubleweights2=oodataframebytimepointreduced$oovecdoubletime %>% purrr::map_dbl(oogetdoublefromdoublesurvivalminusAtt2) %>% purrr::map_dbl(oofunctionweightasafunctionofstminus2)

  #oodataframebytimepointreducedwitht2survival=base::merge(oodataframebytimepointreduced,oodataframeworkingt2,by.x="oovecdoubletime",by.y="time")

  #oovecdoubleweights2=oodataframebytimepointreducedwitht2survival$oovecdoublesurvivalminusATt2 %>% purrr::map_dbl(oofunctionweightasafunctionofstminus2)

  #if(base::length(oovecdoubleweights) != base::length(oovecdoubleweights2) )
  #{
  #  base::warning(base::paste0("oovecdoubleweights length ",base::length(oovecdoubleweights)," and oovecdoubleweights2 length ",base::length(oovecdoubleweights2)," and oovecdoubleweights: ",base::toString(base::round(oovecdoubleweights,digits=2))," and oovecdoubleweights2: ",base::toString(base::round(oovecdoubleweights2,digits = 2) ) ) )
  #}

  #### The above part is new. ####

  #see "Lee 2007 On the versatility of the combination of weighted log-rank statistics.pdf"
  oodoublecovarianceestimate=base::sum( oovecdoubleweights * oovecdoubleweights2 * (oodataframebytimepointreduced$oovecinty1*oodataframebytimepointreduced$oovecinty2)/(oodataframebytimepointreduced$oovecinty1 + oodataframebytimepointreduced$oovecinty2)*(1 - (oodataframebytimepointreduced$oovecintdn1 + oodataframebytimepointreduced$oovecintdn2 - 1)/(oodataframebytimepointreduced$oovecinty1 + oodataframebytimepointreduced$oovecinty2 - 1) ) * (oodataframebytimepointreduced$oovecintdn1 + oodataframebytimepointreduced$oovecintdn2)/(oodataframebytimepointreduced$oovecinty1 + oodataframebytimepointreduced$oovecinty2) )

  oodoublecovarianceestimate

}



oogetdataframemsdataprocessedvRECODED=function(oodataframe,oostringorsymbolid="id",oostringorsymboltreated="treated",oostringorsymbolAtime="Atime",oostringorsymbolBtime="Btime",oostringorsymbolBobserved="Bobserved",oostringorsymbolCtime="Ctime",oostringorsymbolCobserved="Cobserved")
{
  ooexpressionid=rlang::ensym(oostringorsymbolid)
  ooexpressiontreated=rlang::ensym(oostringorsymboltreated)
  ooexpressionAtime=rlang::ensym(oostringorsymbolAtime)
  ooexpressionBtime=rlang::ensym(oostringorsymbolBtime)
  ooexpressionBobserved=rlang::ensym(oostringorsymbolBobserved)
  ooexpressionCtime=rlang::ensym(oostringorsymbolCtime)
  ooexpressionCobserved=rlang::ensym(oostringorsymbolCobserved)

  oomatrixinttmat=base::matrix(NA_integer_,nrow=3L,ncol=3L)
  oomatrixinttmat[[1L,2L]]=1L; oomatrixinttmat[[1L,3L]]=2L;
  base::dimnames(oomatrixinttmat)=base::list(from=c("A","B","C"),to=c("A","B","C") )

  oodataframesupport=oodataframe %>% dplyr::select(-c(!!ooexpressionAtime,!!ooexpressionBtime,!!ooexpressionBobserved,!!ooexpressionCtime,!!ooexpressionCobserved) ) #auxiliary covariates, since mstate::msprep keep= argument seems to be broken
  oodataframemsdata=mstate::msprep(time=c(rlang::as_string(ooexpressionAtime),rlang::as_string(ooexpressionBtime),rlang::as_string(ooexpressionCtime) ),
                                   status=c(NA,rlang::as_string(ooexpressionBobserved),rlang::as_string(ooexpressionCobserved) ),
                                   data=oodataframe,
                                   trans=oomatrixinttmat,
                                   id=rlang::as_string(ooexpressionid),
                                   start=base::list(state=base::rep(1L,base::nrow(oodataframe) ),time=oodataframe[[rlang::as_string(ooexpressionAtime)]] )
  )
  oodataframemsdata=base::merge(oodataframemsdata,oodataframesupport,by=rlang::as_string(ooexpressionid) )
  oodataframemsdata
}




#'
#' @importFrom mvtnorm pmvnorm
#' @importFrom mvtnorm GenzBretz
#'

maxcombo.pvalue <- function(treated, time, status){
  N <- length(treated)
  studyentry <- rep(0,N)
  id <- seq(1:N)
  arm <- 1-treated
  statusC <- 1 - status

  df <- data.frame(id = id, treated = arm, Atime = studyentry, Btime = time, Bobserved = status, Ctime = time, Cobserved = statusC)


  ## MaxCombo package with weights
  oolistfunctionweightasafunctionofstminus=list(
    function(stminus){ base::return(1.0) }, # (0,0)
    function(stminus){ base::return(1.0 - stminus) }, # (0,1)
    function(stminus){ base::return(stminus) } # (1,0)
  )


  # oogetdoublemaxcombotestpvalue - Function
  ooexpressionid="id"
  ooexpressiontreated="treated"
  ooexpressionAtime="Atime"
  ooexpressionBtime="Btime"
  ooexpressionBobserved="Bobserved"
  ooexpressionCtime="Ctime"
  ooexpressionCobserved="Cobserved"

  oodoublemaxcomboteststatistic = oogetdoublemaxcomboteststatistic(oodataframe = df,oolistfunctionweightasafunctionofstminus = oolistfunctionweightasafunctionofstminus,oostringorsymbolid = !!ooexpressionid,oostringorsymboltreated = !!ooexpressiontreated,oostringorsymbolAtime = !!ooexpressionAtime,oostringorsymbolBtime = !!ooexpressionBtime,oostringorsymbolBobserved = !!ooexpressionBobserved,oostringorsymbolCtime = !!ooexpressionCtime,oostringorsymbolCobserved = !!ooexpressionCobserved)

  oointp = base::length(oolistfunctionweightasafunctionofstminus)
  oomatrixdoublecovariancematrix = base::matrix(NA,nrow = oointp,ncol = oointp)

  for(oointj in 1L:oointp)
  {
    for(oointk in 1L:oointj)
    {
      oomatrixdoublecovariancematrix[[oointj,oointk]]=oogetresultsbasiccovariance3vRECODED(oodataframet1 = df,oodataframet2 = df,oofunctionweightasafunctionofstminus1 = oolistfunctionweightasafunctionofstminus[[oointj]],oofunctionweightasafunctionofstminus2 = oolistfunctionweightasafunctionofstminus[[oointk]],oostringorsymbolid = !!ooexpressionid,oostringorsymboltreated = !!ooexpressiontreated,oostringorsymbolAtime = !!ooexpressionAtime,oostringorsymbolBtime = !!ooexpressionBtime,oostringorsymbolBobserved = !!ooexpressionBobserved,oostringorsymbolCtime = !!ooexpressionCtime,oostringorsymbolCobserved = !!ooexpressionCobserved)
      oomatrixdoublecovariancematrix[[oointk,oointj]]=oomatrixdoublecovariancematrix[[oointj,oointk]]
    }
  }

  oomatrixdoublecorrelationmatrix = oomatrixdoublecovariancematrix %>% stats::cov2cor()

  lower <- -abs(rep(oodoublemaxcomboteststatistic,times=oointp))
  upper <- abs(rep(oodoublemaxcomboteststatistic,times=oointp))
  mean <- rep(0, times = oointp)
  two_tailed_p_value = 1 - pmvnorm(mean=mean,sigma=oomatrixdoublecorrelationmatrix,lower=lower,upper=upper, algorithm= GenzBretz(maxpts=50000,abseps=0.00001))
  return(two_tailed_p_value)
}
