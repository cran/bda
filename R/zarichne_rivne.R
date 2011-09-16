
getbdata <- function(x,type='head',gender="male",location="rc",filter=FALSE){
  name <- deparse(substitute(x))
  if(filter){
    sele1 = x$Gestation >= 38
    sele2 = x$Plurality == 1
    sele3 = x$Live.Still == 'liveborn'
    x = x[sele1&sele2&sele3,]
    ###  The following filters will be applied individually, in
    ###  regression analysis, we need consider apply these filter as
    ###  well.
    
    ##    sele4 = x$Head > 0
    ##    sele5 = x$Mother.s.age>0 & x$Mother.s.age<60
    ##    sele6 = x$Weight > 0
    ##    x = x[sele1&sele2&sele3&sele4&sele5&sele6,]
  }
  location = match.arg(tolower(location),c("r","u","zz","bz","sz","rc"))
  gender = match.arg(tolower(gender),c("male","female"))
  type = match.arg(tolower(type),
    c("length","weight","head","ofc","chest","mother.s.age","ga","gestation"))
  sele1 = tolower(x$Sex) == gender
  sele2 = x$region == location
  sele3 = x$type == location
  sele = sele1 & (sele2 | sele3)
  if(type=="length"){
    y = x$Length[sele]
    just = 'center'
    scale = 1
  }else if(type=="weight"){
    y = x$Weight[sele]
    just = 'center'
    scale = 100
  }else if(type=="head"||type=="ofc"){
    y = x$Head[sele]
    just = 'center'
    scale = 1
  }else if(type=="chest"){y
    y = x$Chest[sele]
    just = 'center'
    scale = 1
  }else if(type=="mother.s.age"){
    y = x$Mother.s.age[sele]
    just = 'center'
    scale = 1
  }else if(type=="ga"||type=="gestation"){
    y = x$Gestation[sele]
    just = 'left'
    scale = 1
  }
  x0 = .discretize(y,na.rm=TRUE,just=just)
  structure(list(x0,
                 name = paste(name,"(",type,",scale=",scale,")",sep='')
                 ), class='bdata')
}

##  out = switch(method,
##    em = binem(y,m=m,mu=mu,scal=scale,just=just),
##    smkde = smkde(y,scale=scale,just=just),
##    histospline = histospline(y,scale=scale,just=just)
##  )
##  out

