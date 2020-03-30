
mediation.test <- function(mv,iv,dv){
    ## mx: mediation variable
    ## iv:  indep. variable
    ## dv: dep. var.
    nm = length(mv); ni = length(iv); nd = length(dv);
    if(nm!=ni | nm!=nd | ni!=nd) stop("Variables have different lengths.")

    sele1 <- is.na(mv)
    sele2 <- is.na(iv)
    sele3 <- is.na(dv)
    sele <- sele1 | sele2 | sele3
    nmi <- sum(sele)
    if(nmi > 0){
        warning("missing value(s) removed")
    }
    
    if(nm==nmi){
        stop("all records contain missing values")
    }else if(nm - nmi < 5){
        stop("too few valid records to test")
    }else if(nm - nmi < 10){
        warning("the number of records might be too small")
    }
    
    tmp = summary(lm(mv~iv));
    a = tmp$coef[2,1];sa=tmp$coef[2,2];
    tmp = summary(lm(dv~mv+iv));
    b = tmp$coef[2,1];sb=tmp$coef[2,2];
    tmp1 = b^2*sa^2+a^2*sb^2
    tmp2 = sa^2*sb^2
    zsob = a*b/sqrt(tmp1);
    psob = pnorm(-abs(zsob))*2;
    zaro = a*b/sqrt(tmp1+tmp2);
    paro = pnorm(-abs(zaro))*2;
    if(tmp1>tmp2){
        zgm = a*b/sqrt(tmp1-tmp2)
        pgm = pnorm(-abs(zgm))*2;
    }else{
        zgm=NA
        pgm=NA;
    }
    p.value=c(psob,paro,pgm)
    z.value = c(zsob,zaro,zgm)
    out = data.frame(rbind(z.value,p.value));
    names(out)=c("Sobel","Aroian","Goodman")
    out
}
