`PKmodelFinder` <-
function()
{  options(guiToolkit="RGtk2")
   require(gWidgets)
   require(tcltk)
   require(tkrplot)
   require(RGtk2)
   require(gWidgetsRGtk2)
   require(cairoDevice)
   require(numDeriv)

   
   ### One Comp First Order Absorption

   Comp1.FO.L<-function(dose,Tinf,time,mu)
   {  logV<-mu[1]  ; V<-exp(logV)
      logKa<-mu[2] ; Ka<-exp(logKa)
      logKe<-mu[3] ; Ke<-exp(logKe)
      dose*Ka*(exp(-Ke*time)-exp(-Ka*time))/(V*(Ka-Ke))
   }

   I.Comp1.FO.L<-function(dose,Tinf,time,logV,logKa,logKe)
   {  dose*exp(logKa)*(exp(-exp(logKe)*time)-
      exp(-exp(logKa)*time))/(exp(logV)*(exp(logKa)-exp(logKe)))
   }

   I.Comp1.FO.L.optim<-function(x)
   {  logV<-x[1]
      logKa<-x[2]
      logKe<-x[3]
      ID.list<-unique(tot.data$subject)
      I.value<-0
      for(i in ID.list)
      {  id.s<-which(tot.data$subject==i)
         temp1<- tot.data$conc[id.s]-I.Comp1.FO.L(tot.data$dose[id.s],tot.data$Tinf[id.s],
                                              tot.data$time[id.s],logV,logKa,logKe)  
         I.value<-I.value+sum(temp1*temp1)
       }
       I.value
    }

    One.FirstOrder.linear<-function(f,f.o,a,b,n.theta,model)
    {   dataMake(a,b)
        setModel(n.theta,model)
        setInitial(f,f.o,model)
    }

    One.FirstOrder.linear1<-function(f,f.i,param.list) 
    {   Simulation(f,m,K)
        NDfunction<-"numericDeriv(quote(loglike.ind(id.log,T.data,f,
                              logV,logKa,logKe)),c('logV','logKa','logKe'), myenv)"
        FisherInfo(f.i,param.list,NDfunction)
        Pred(f)
        return(tot.data)
    }
    
   ### One Comp IV Bolus  	
     Comp1.IVbolus.L<-function(dose,Tinf,time,mu)
   {  logV<-mu[1]
      logKe<-mu[2]
      dose*(exp(-exp(logKe)*time))/exp(logV)
   }

   I.Comp1.IVbolus.L<-function(dose,Tinf,time,logV,logKe)
   {  dose*(exp(-exp(logKe)*time))/exp(logV)
   }

   I.Comp1.IVbolus.L.optim<-function(x)
   {  logV<-x[1]
      logKe<-x[2]
      ID.list<-unique(tot.data$subject)
      I.value<-0
      for(i in ID.list)
      {  id.s<-which(tot.data$subject==i)
          temp1<- tot.data$conc[id.s]-I.Comp1.IVbolus.L(tot.data$dose[id.s],tot.data$Tinf[id.s],
                                      tot.data$time[id.s],logV,logKe)  
          I.value<-I.value+sum(temp1*temp1)
      }
      I.value
   }

   One.IVbolus.linear<-function(f,f.o,a,b,n.theta,model)
   {  dataMake(a,b)
      setModel(n.theta,model)
      setInitial(f,f.o,model)
   }
   
   One.IVbolus.linear1<-function(f,f.i,param.list)
   {  Simulation(f,m,K)
      NDfunction<-"numericDeriv(quote(loglike.ind(id.log,T.data,f,
                logV,logKe)),c('logV','logKe'), myenv)"
      FisherInfo(f.i,param.list,NDfunction)
      Pred(f)
      return(tot.data)
   }
   ### One Comp IV infusion
      Comp1.IVinf.L<-function(dose,TTinf,time,mu)
   {  logV<-mu[1]
      logKe<-mu[2]
      Tinf<-TTinf[1]
      ifelse(time<=Tinf,
	         dose*(1-exp(-exp(logKe)*time))/(Tinf*exp(logKe)*exp(logV)),
	         dose*(1-exp(-exp(logKe)*Tinf))*exp(-exp(logKe)*(time-Tinf))/(Tinf*exp(logKe)*exp(logV)))  
   }

   I.Comp1.IVinf.L<-function(dose,TTinf,time,logV,logKe)
   {  Tinf<-TTinf[1]
      ifelse(time<=Tinf,
	         dose*(1-exp(-exp(logKe)*time))/(Tinf*exp(logKe)*exp(logV)),
	         dose*(1-exp(-exp(logKe)*Tinf))*exp(-exp(logKe)*(time-Tinf))/(Tinf*exp(logKe)*exp(logV))) 
   }

   I.Comp1.IVinf.L.optim<-function(x)
   {  logV<-x[1]
      logKe<-x[2]
      ID.list<-unique(tot.data$subject)
      I.value<-0
      for(i in ID.list)
      {  id.s<-which(tot.data$subject==i)
         temp1<- tot.data$conc[id.s]-I.Comp1.IVinf.L(tot.data$dose[id.s],tot.data$Tinf[id.s],
                                      tot.data$time[id.s],logV,logKe)  
         I.value<-I.value+sum(temp1*temp1)
      }
      I.value
   }

   One.IVinf.linear<-function(f,f.o,a,b,n.theta,model)
   {  dataMake(a,b)
      setModel(n.theta,model)
      setInitial(f,f.o,model)
   }

   One.IVinf.linear1<-function(f,f.i,param.list)
   {  Simulation(f,m,K)
      NDfunction<-"numericDeriv(quote(loglike.ind(id.log,T.data,f,
                      logV,logKe)),c('logV','logKe'), myenv)"
      FisherInfo(f.i,param.list,NDfunction)
      Pred(f)
      return(tot.data)
   }
   ### One Comp zero order absorption
      Comp1.ZO.L<-function(dose,Tinf,time,mu)
   {  logV<-mu[1]  ; V<-exp(logV)
      logKe<-mu[2] ; Ke<-exp(logKe)
      ifelse(time<=Tk0,
	         (dose/(Tk0*Ke*V))*(1-exp(-Ke*time)),
           (dose/(Tk0*Ke*V))*(1-exp(-Ke*Tk0))*exp(-Ke*(time-Tk0)) )
   }

   I.Comp1.ZO.L<-function(dose,Tinf,time,logV,logKe)
   {  V<-exp(logV)
      Ke<-exp(logKe)
      ifelse(time<=Tk0,
	         (dose/(Tk0*Ke*V))*(1-exp(-Ke*time)),
	         (dose/(Tk0*Ke*V))*(1-exp(-Ke*Tk0))*exp(-Ke*(time-Tk0)) )
   }

   I.Comp1.ZO.L.optim<-function(x)
   {  logV<-x[1]
      logKe<-x[2]
      ID.list<-unique(tot.data$subject)
      I.value<-0
      for(i in ID.list)
      {  id.s<-which(tot.data$subject==i)
         temp1<- tot.data$conc[id.s]-I.Comp1.ZO.L(tot.data$dose[id.s],tot.data$Tinf[id.s],
                                      tot.data$time[id.s],logV,logKe)  
         I.value<-I.value+sum(temp1*temp1)
      }
      I.value
   }

   One.ZeroOrder.linear<-function(f,f.o,a,b,n.theta,model)
   {  dataMake(a,b)
      setModel(n.theta,model)
      setInitial(f,f.o,model)
   }

   One.ZeroOrder.linear1<-function(f,f.i,param.list)
   {  Simulation(f,m,K)
      NDfunction<-"numericDeriv(quote(loglike.ind(id.log,T.data,f,
                         logV,logKe)),c('logV','logKe'), myenv)"
      FisherInfo(f.i,param.list,NDfunction)
      Pred(f)
      return(tot.data)
  }
   ### two comp first order absorption
      Comp2.FO.L<-function(dose,Tinf,time,mu)
   {  logV<-mu[1] ; V<-exp(logV)
      logKa<-mu[2]; Ka<-exp(logKa)
      logKe<-mu[3]; Ke<-exp(logKe)
      logK12<-mu[4]; K12<-exp(logK12)
      logK21<-mu[5]; K21<-exp(logK21)

      beta<-0.5*(K12+K21+Ke-sqrt((K12+K21+Ke)^2-4*K21*Ke))
      alpha<-K21*Ke/beta
      A<-Ka*(K21-alpha)/(V*(Ka-alpha)*(beta-alpha))
      B<-Ka*(K21-beta)/(V*(Ka-beta)*(alpha-beta))  
      dose*(A*exp(-alpha*time)+B*exp(-beta*time)-(A+B)*exp(-Ka*time))
   }

   I.Comp2.FO.L<-function(dose,Tinf,time,logV,logKa,logKe,logK12,logK21)
   {  V<-exp(logV)
      Ka<-exp(logKa)
      K12<-exp(logK12)
      K21<-exp(logK21)
      Ke<-exp(logKe)
      beta<-0.5*(K12+K21+Ke-sqrt((K12+K21+Ke)^2-4*K21*Ke))
      alpha<-K21*Ke/beta
      A<-Ka*(K21-alpha)/(V*(Ka-alpha)*(beta-alpha))
      B<-Ka*(K21-beta)/(V*(Ka-beta)*(alpha-beta))  
      dose*(A*exp(-alpha*time)+B*exp(-beta*time)-(A+B)*exp(-Ka*time))
   }

   I.Comp2.FO.L.optim<-function(x)
   {  logV<-x[1] ; 
      logKa<-x[2]; 
      logK12<-x[3]; 
      logK21<-x[4]; 
      logKe<-x[5]; 
      ID.list<-unique(tot.data$subject)
      I.value<-0
      for(i in ID.list)
      {  id.s<-which(tot.data$subject==i)
         temp1<- tot.data$conc[id.s]-I.Comp2.FO.L(tot.data$dose[id.s],tot.data$Tinf[id.s],
                                      tot.data$time[id.s],logV,logKa,logKe,logK12,logK21)  
         I.value<-I.value+sum(temp1*temp1)
      }
      I.value
   }
    
   Two.FirstOrder.linear<-function(f,f.o,a,b,n.theta,model)
   {  dataMake(a,b)
      setModel(n.theta,model)
      setInitial(f,f.o,model)
   }
    
   Two.FirstOrder.linear1<-function(f,f.i,param.list)
   {  Simulation(f,m,K)
      NDfunction<-"numericDeriv(quote(loglike.ind(id.log,T.data,f,logV,logKa,logKe,logK12,logK21)),
                                c('logV','logKa','logKe','logK12','logK21'),Keep.env)"
      FisherInfo(f.i,param.list,NDfunction)
      Pred(f)
      return(tot.data)
   }
   ### two Comp IV Bolus  
      Comp2.IVbolus.L<-function(dose,Tinf,time,mu)
   {  logV<-mu[1]; V<-exp(logV)
      logKe<-mu[2]; Ke<-exp(logKe)
      logK12<-mu[3]; K12<-exp(logK12)
      logK21<-mu[4]; K21<-exp(logK21)
      beta<-0.5*(K12+K21+Ke-sqrt((K12+K21+Ke)^2-4*K21*Ke))
      alpha<-K21*Ke/beta
      A<-(alpha-K21)/V*(alpha-beta)
      B<-(beta-K21)/V*(beta-alpha)
      dose*(A*exp(-alpha*time)+B*exp(-beta*time))
   }

   I.Comp2.IVbolus.L<-function(dose,Tinf,time,logV,logKe,logK12,logK21)
   {  V<-exp(logV)
      Ke<-exp(logKe)
      K12<-exp(logK12)
      K21<-exp(logK21)
      beta<-0.5*(K12+K21+Ke-sqrt((K12+K21+Ke)^2-4*K21*Ke))
      alpha<-K21*Ke/beta
      A<-(alpha-K21)/V*(alpha-beta)
      B<-(beta-K21)/V*(beta-alpha)
      dose*(A*exp(-alpha*time)+B*exp(-beta*time))
   }

   I.Comp2.IVbolus.L.optim<-function(x)
   {  logV<-x[1]; 
      logKe<-x[2]; 
      logK12<-x[3]; 
      logK21<-x[4]; 
      ID.list<-unique(tot.data$subject)
      I.value<-0
      for(i in ID.list)
      {   id.s<-which(tot.data$subject==i)
          temp1<- tot.data$conc[id.s]-I.Comp2.IVbolus.L(tot.data$dose[id.s],tot.data$Tinf[id.s],
                                      tot.data$time[id.s],logV,logKe,logK12,logK21)  
          I.value<-I.value+sum(temp1*temp1)
      }
      I.value
   }

   Two.IVbolus.linear<-function(f,f.o,a,b,n.theta,model)
   {  dataMake(a,b)
      setModel(n.theta,model)
      setInitial(f,f.o,model)
   }

   Two.IVbolus.linear1<-function(f,f.i,param.list)
   {  Simulation(f,m,K)
      NDfunction<-"numericDeriv(quote(loglike.ind(id.log,T.data,f,logV,logKe,logK12,logK21)),
                      c('logV','logKe','logK12','logK21'), myenv)"
      FisherInfo(f.i,param.list,NDfunction)
      Pred(f)
   }  
   ### two Comp IV infusion
      Comp2.IVinf.L<-function(dose,TTinf,time,mu)
   {  Tinf<-TTinf[1]
      logV<-mu[1]; V<-exp(logV)
      logKe<-mu[2]; Ke<-exp(logKe)
      logK12<-mu[3]; K12<-exp(logK12)
      logK21<-mu[4]; K21<-exp(logK21)
      beta<-0.5*(K12+K21+Ke-sqrt((K12+K21+Ke)^2-4*K21*Ke))
      alpha<-K21*Ke/beta
      A<-(alpha-K21)/V*(alpha-beta)
      B<-(beta-K21)/V*(beta-alpha)
      ifelse(time<=Tinf,
      	(dose/Tinf)*((A/alpha)*(1-exp(-alpha*time))+(B/beta)*(1-exp(-beta*time))),
      	(dose/Tinf)*((A/alpha)*(1-exp(-alpha*Tinf))*exp(-alpha*(time-Tinf))+
                                 (B/beta)*(1-exp(-beta*Tinf))*exp(-beta*(time-Tinf)))) 
   }

   I.Comp2.IVinf.L<-function(dose,TTinf,time,logV,logKe,logK12,logK21)
   {  Tinf<-TTinf[1]
      V<-exp(logV)
      Ke<-exp(logKe)
      K12<-exp(logK12)
      K21<-exp(logK21)
      beta<-0.5*(K12+K21+Ke-sqrt((K12+K21+Ke)^2-4*K21*Ke))
      alpha<-K21*Ke/beta
      A<-(alpha-K21)/V*(alpha-beta)
      B<-(beta-K21)/V*(beta-alpha)
      ifelse(time<=Tinf,
      	(dose/Tinf)*((A/alpha)*(1-exp(-alpha*time))+(B/beta)*(1-exp(-beta*time))),
      	(dose/Tinf)*((A/alpha)*(1-exp(-alpha*Tinf))*exp(-alpha*(time-Tinf))+
                                 (B/beta)*(1-exp(-beta*Tinf))*exp(-beta*(time-Tinf))))
   }

   I.Comp2.IVinf.L.optim<-function(x)
   {  logV<-x[1]; 
      logKe<-x[2]; 
      logK12<-x[3];
      logK21<-x[4];
      ID.list<-unique(tot.data$subject)
      I.value<-0
      for(i in ID.list)
      {  id.s<-which(tot.data$subject==i)
         temp1<- tot.data$conc[id.s]-I.Comp2.IVinf.L(tot.data$dose[id.s],tot.data$Tinf[id.s],
                                       tot.data$time[id.s],logV,logKe,logK12,logK21)  
         I.value<-I.value+sum(temp1*temp1)
      }
      I.value
   }

   Two.IVinf.linear<-function(f,f.o,a,b,n.theta,model)
   {   dataMake(a,b)
       setModel(n.theta,model)
       setInitial(f,f.o,model)
   }

   Two.IVinf.linear1<-function(f,f.i,param.list)
   {  Simulation(f,m,K)
      NDfunction<-"numericDeriv(quote(loglike.ind(id.log,T.data,f,logV,logKe,logK12,logK21)),c('logV','logKe','logK12','logK21'), myenv)"
      FisherInfo(f.i,param.list,NDfunction)
      Pred(f)
      return(tot.data)
   }
   ### two Comp zero order absorption 
      Comp2.ZO.L<-function(dose,Tinf,time,mu)
   {  logV<-mu[1]  ; V<-exp(logV)
      logKe<-mu[2] ; Ke<-exp(logKe)
      logK12<-mu[3] ; K12<-exp(logK12)
      logK21<-mu[4] ; K21<-exp(logK21)
      beta<-0.5*(K12+K21+Ke-sqrt((K12+K21+Ke)^2-4*K21*Ke))
      alpha<-K21*Ke/beta
      A<-(alpha-K21)/(V*(alpha-beta))
      B<-(beta-K21)/(V*(beta-alpha))

      ifelse(time<=Tk0,
 	      (dose/Tk0)*( (A/alpha)*(1-exp(-alpha*time))+(B/beta)*(1-exp(-beta*time)) ),
	       (dose/Tk0)*( (A/alpha)*(1-exp(-alpha*Tk0))*exp(-alpha*(time-Tk0))+(B/beta)*(1-exp(-beta*Tk0))*exp(-alpha*(time-Tk0)) ) ) 
   }

   I.Comp2.ZO.L<-function(dose,Tinf,time,logV,logKe,logK12,logK21)
   {  V<-exp(logV)
      Ke<-exp(logKe)
      K12<-exp(logK12)
      K21<-exp(logK21)

      beta<-0.5*(K12+K21+Ke-sqrt((K12+K21+Ke)^2-4*K21*Ke))
      alpha<-K21*Ke/beta
      A<-(alpha-K21)/(V*(alpha-beta))
      B<-(beta-K21)/(V*(beta-alpha))

      ifelse(time<=Tk0,
 	        (dose/Tk0)*( (A/alpha)*(1-exp(-alpha*time))+(B/beta)*(1-exp(-beta*time)) ),
	        (dose/Tk0)*( (A/alpha)*(1-exp(-alpha*Tk0))*exp(-alpha*(time-Tk0))+(B/beta)*(1-exp(-beta*Tk0))*exp(-alpha*(time-Tk0)) ) )
   }

   I.Comp2.ZO.L.optim<-function(x)
   {  logV<-x[1]
      logKe<-x[2]
      logK12<-x[3]
      logK21<-x[4]
      ID.list<-unique(tot.data$subject)
      I.value<-0
      for(i in ID.list)
      {  id.s<-which(tot.data$subject==i)
         temp1<- tot.data$conc[id.s]-I.Comp2.ZO.L(tot.data$dose[id.s],tot.data$Tinf[id.s],
                                      tot.data$time[id.s],logV,logKe,logK12,logK21)  
        I.value<-I.value+sum(temp1*temp1)
      }
      I.value
   }

   Two.ZeroOrder.linear<-function(f,f.o,a,b,n.theta,model)
   {  dataMake(a,b)
      setModel(n.theta,model)
      setInitial(f,f.o,model)
   }

   Two.ZeroOrder.linear1<-function(f,f.i,param.list)
   {  Simulation(f,m,K)
      NDfunction<-"numericDeriv(quote(loglike.ind(id.log,T.data,f,logV,logKe,logK12,logK21)),c('logV','logKe','logK12','logK21'), myenv)"
      FisherInfo(f.i,param.list,NDfunction)
      Pred(f)
      return(tot.data)
   }
   ### three comp first order absorption
   Comp3.FO.L<-function(dose,Tinf,time,mu)
   {  logV<-mu[1] ; V<-exp(logV)
      logKa<-mu[2]; Ka<-exp(logKa)
      logKe<-mu[3]; Ke<-exp(logKe)
      logK12<-mu[4]; K12<-exp(logK12)
      logK21<-mu[5]; K21<-exp(logK21)
      logK13<-mu[6]; K13<-exp(logK13)
      logK31<-mu[7]; K31<-exp(logK31)

      a0<-Ke*K21*K31;a1<-Ke*K31+K21*K31+K21*K13+Ke*K21+K31*K12;a2<-Ke+K12+K13+K21+K31
      p<-a1-(a2^2)/3;q<-((2*a2^3)/27)-(a1*a2/3)+a0
      r1<-sqrt(-(p^3/27));r2<-2*r1^(1/3);PI<-acos(-q/(2*r1))/3
      alpha<--(cos(PI)*r2-(a2/3))
      beta<--(cos(PI+(2*pi/3))*r2-(a2/3))
      gamma<--(cos(PI+(4*pi/3))*r2-(a2/3))
      A<-(Ka*(K21-alpha)*(K31-alpha))/(V*(Ka-alpha)*(alpha-beta)*(alpha-gamma))
      B<-(Ka*(K21-beta)*(K31-beta))/(V*(Ka-beta)*(beta-alpha)*(beta-gamma))
      C<-(Ka*(K21-gamma)*(K31-gamma))/(V*(Ka-gamma)*(gamma-beta)*(gamma-alpha))
      dose*(A*exp(-alpha*T)+B*exp(-beta*time)+C*exp(-gamma*time)-(A+B+C)*exp(-Ka*time))
  }

  I.Comp3.FO.L<-function(dose,Tinf,time,logV,logKa,logKe,logK12,logK21,logK13,logK31)
  {   V<-exp(logV)
      Ka<-exp(logKa)
      K12<-exp(logK12)
      K21<-exp(logK21)
      K13<-exp(logK13)
      K31<-exp(logK31)
      Ke<-exp(logKe)
      a0<-Ke*K21*K31;a1<-Ke*K31+K21*K31+K21*K13+Ke*K21+K31*K12;a2<-Ke+K12+K13+K21+K31
      p<-a1-(a2^2)/3;q<-((2*a2^3)/27)-(a1*a2/3)+a0
      r1<-sqrt(-(p^3/27));r2<-2*r1^(1/3);PI<-acos(-q/(2*r1))/3
      alpha<--(cos(PI)*r2-(a2/3))
      beta<--(cos(PI+(2*pi/3))*r2-(a2/3))
      gamma<--(cos(PI+(4*pi/3))*r2-(a2/3))
      A<-(Ka*(K21-alpha)*(K31-alpha))/(V*(Ka-alpha)*(alpha-beta)*(alpha-gamma))
      B<-(Ka*(K21-beta)*(K31-beta))/(V*(Ka-beta)*(beta-alpha)*(beta-gamma))
      C<-(Ka*(K21-gamma)*(K31-gamma))/(V*(Ka-gamma)*(gamma-beta)*(gamma-alpha))
      dose*(A*exp(-alpha*T)+B*exp(-beta*time)+C*exp(-gamma*time)-(A+B+C)*exp(-Ka*time))
  }

  I.Comp3.FO.L.optim<-function(mu)
  {  logV<-mu[1] ; 
     logKa<-mu[2]; 
     logKe<-mu[3]; 
     logK12<-mu[4]; 
     logK21<-mu[5]; 
     logK13<-mu[6]; 
     logK31<-mu[7]; 
     ID.list<-unique(tot.data$subject)
     I.value<-0
     for(i in ID.list)
     {  id.s<-which(tot.data$subject==i)
        temp1<- tot.data$conc[id.s]-I.Comp3.FO.L(tot.data$dose[id.s],tot.data$Tinf[id.s],
                                tot.data$time[id.s],logV,logKa,logKe,logK12,logK21,logK13,logK31)  
        I.value<-I.value+sum(temp1*temp1)
     }
     I.value
  }
    
  Three.FirstOrder.linear<-function(f,f.o,a,b,n.theta,model)
  {  dataMake(a,b)
     setModel(n.theta,model)
     setInitial(f,f.o,model)
  }
    
   Three.FirstOrder.linear1<-function(f,f.i,param.list)
  {  Simulation(f,m,K)
     NDfunction<-"numericDeriv(quote(loglike.ind(id.log,T.data,f,logV,logKa,logKe,logK12,logK21,logK13,logK31)),
                                c('logV','logKa','logKe','logK12','logK21','logK13','logK31'),Keep.env)"
     FisherInfo(f.i,param.list,NDfunction)
     Pred(f)
     return(tot.data)
  }

   
   ### three Comp IV Bolus  
   Comp3.IVbolus.L<-function(dose,Tinf,time,mu)
   {  logV<-mu[1] ; V<-exp(logV)
      logKe<-mu[2]; Ke<-exp(logKe)
      logK12<-mu[3]; K12<-exp(logK12)
      logK21<-mu[4]; K21<-exp(logK21)
      logK13<-mu[5]; K13<-exp(logK13)
      logK31<-mu[6]; K31<-exp(logK31)
      a0<-Ke*K21*K31;a1<-Ke*K31+K21*K31+K21*K13+Ke*K21+K31*K12;a2<-Ke+K12+K13+K21+K31
      p<-a1-(a2^2)/3;q<-((2*a2^3)/27)-(a1*a2/3)+a0
      r1<-sqrt(-(p^3/27));r2<-2*r1^(1/3);PI<-acos(-q/(2*r1))/3
      alpha<--(cos(PI)*r2-(a2/3))
      beta<--(cos(PI+(2*pi/3))*r2-(a2/3))
      gamma<--(cos(PI+(4*pi/3))*r2-(a2/3))
      A<-((K21-alpha)*(K31-alpha))/(V*(alpha-beta)*(alpha-gamma))
      B<-((K21-beta)*(K31-beta))/(V*(beta-alpha)*(beta-gamma))
      C<-((K21-gamma)*(K31-gamma))/(V*(gamma-beta)*(gamma-alpha))
      dose*(A*exp(-alpha*time)+B*exp(-beta*time)+C*exp(-gamma*time))
   }

   I.Comp3.IVbolus.L<-function(dose,Tinf,time,logV,logKe,logK12,logK21,logK13,logK31)
   {  V<-exp(logV)
      Ke<-exp(logKe)
      K12<-exp(logK12)
      K21<-exp(logK21)
      K13<-exp(logK13)
      K31<-exp(logK31)
      a0<-Ke*K21*K31;a1<-Ke*K31+K21*K31+K21*K13+Ke*K21+K31*K12;a2<-Ke+K12+K13+K21+K31
      p<-a1-(a2^2)/3;q<-((2*a2^3)/27)-(a1*a2/3)+a0
      r1<-sqrt(-(p^3/27));r2<-2*r1^(1/3);PI<-acos(-q/(2*r1))/3
      alpha<--(cos(PI)*r2-(a2/3))
      beta<--(cos(PI+(2*pi/3))*r2-(a2/3))
      gamma<--(cos(PI+(4*pi/3))*r2-(a2/3))
      A<-((K21-alpha)*(K31-alpha))/(V*(alpha-beta)*(alpha-gamma))
      B<-((K21-beta)*(K31-beta))/(V*(beta-alpha)*(beta-gamma))
      C<-((K21-gamma)*(K31-gamma))/(V*(gamma-beta)*(gamma-alpha))
      dose*(A*exp(-alpha*time)+B*exp(-beta*time)+C*exp(-gamma*time))
   }

   I.Comp3.IVbolus.L.optim<-function(x)
   {  logV<-x[1]; 
      logKe<-x[2]; 
      logK12<-x[3]; 
      logK21<-x[4]; 
      logK13<-x[5]; 
      logK31<-x[6]; 
      ID.list<-unique(tot.data$subject)
      I.value<-0
      for(i in ID.list)
      {  id.s<-which(tot.data$subject==i)
         temp1<- tot.data$conc[id.s]-I.Comp3.IVbolus.L(tot.data$dose[id.s],tot.data$Tinf[id.s],
                                       tot.data$time[id.s],logV,logKe,logK12,logK21,logK13,logK31)  
         I.value<-I.value+sum(temp1*temp1)
      }
      I.value
   }
   
   Three.IVbolus.linear<-function(f,f.o,a,b,n.theta,model)
   {  dataMake(a,b)
      setModel(n.theta,model)
      setInitial(f,f.o,model)
   }
   
   Three.IVbolus.linear1<-function(f,f.i,param.list)
   {  Simulation(f,m,K)
      NDfunction<-"numericDeriv(quote(loglike.ind(id.log,T.data,f,logV,logKe,logK12,logK21,logK13,logK31)),
                      c('logV','logKe','logK12','logK21','logK13','logK31'), myenv)"
      FisherInfo(f.i,param.list,NDfunction)
      Pred(f)
      return(tot.data)
   }

     
   ### three Comp IV infusion
   Comp3.IVinf.L<-function(dose,TTinf,time,mu)
   {  Tinf<-TTinf[1]
      logV<-mu[1]; V<-exp(logV)
      logKe<-mu[2]; Ke<-exp(logKe)
      logK12<-mu[3]; K12<-exp(logK12)
      logK21<-mu[4]; K21<-exp(logK21)
      logK13<-mu[5]; K13<-exp(logK13)
      logK31<-mu[6]; K31<-exp(logK31)

      a0<-Ke*K21*K31;a1<-Ke*K31+K21*K31+K21*K13+Ke*K21+K31*K12;a2<-Ke+K12+K13+K21+K31
      p<-a1-(a2^2)/3;q<-((2*a2^3)/27)-(a1*a2/3)+a0
      r1<-sqrt(-(p^3/27));r2<-2*r1^(1/3);PI<-acos(-q/(2*r1))/3
      alpha<--(cos(PI)*r2-(a2/3))
      beta<--(cos(PI+(2*pi/3))*r2-(a2/3))
      gamma<--(cos(PI+(4*pi/3))*r2-(a2/3))
      A<-((K21-alpha)*(K31-alpha))/(V*(alpha-beta)*(alpha-gamma))
      B<-((K21-beta)*(K31-beta))/(V*(beta-alpha)*(beta-gamma))
      C<-((K21-gamma)*(K31-gamma))/(V*(gamma-beta)*(gamma-alpha))
      ifelse(time<=Tinf,
      	(dose/Tinf)*((A/alpha)*(1-exp(-alpha*time))+
                                 (B/beta)*(1-exp(-beta*time))+
                                 (C/gamma)*(1-exp(-gamma*time))),
      	(dose/Tinf)*((A/alpha)*(1-exp(-alpha*Tinf))*exp(-alpha*(time-Tinf))+
                                 (B/beta)*(1-exp(-beta*Tinf))*exp(-beta*(time-Tinf))+
                                 (C/gamma)*(1-exp(-gamma*Tinf))*exp(-gamma*(time-Tinf)))) 
   }

   I.Comp3.IVinf.L<-function(dose,TTinf,time,logV,logKe,logK12,logK21,logK13,logK31)
   {  Tinf<-TTinf[1]
      V<-exp(logV)
      Ke<-exp(logKe)
      K12<-exp(logK12)
      K21<-exp(logK21)
      K13<-exp(logK13)
      K31<-exp(logK31)

      a0<-Ke*K21*K31;a1<-Ke*K31+K21*K31+K21*K13+Ke*K21+K31*K12;a2<-Ke+K12+K13+K21+K31
      p<-a1-(a2^2)/3;q<-((2*a2^3)/27)-(a1*a2/3)+a0
      r1<-sqrt(-(p^3/27));r2<-2*r1^(1/3);PI<-acos(-q/(2*r1))/3
      alpha<--(cos(PI)*r2-(a2/3))
      beta<--(cos(PI+(2*pi/3))*r2-(a2/3))
      gamma<--(cos(PI+(4*pi/3))*r2-(a2/3))
      A<-((K21-alpha)*(K31-alpha))/(V*(alpha-beta)*(alpha-gamma))
      B<-((K21-beta)*(K31-beta))/(V*(beta-alpha)*(beta-gamma))
      C<-((K21-gamma)*(K31-gamma))/(V*(gamma-beta)*(gamma-alpha))
      ifelse(time<=Tinf,
      	(dose/Tinf)*((A/alpha)*(1-exp(-alpha*time))+
                                 (B/beta)*(1-exp(-beta*time))+
                                 (C/gamma)*(1-exp(-gamma*time))),
      	(dose/Tinf)*((A/alpha)*(1-exp(-alpha*Tinf))*exp(-alpha*(time-Tinf))+
                                 (B/beta)*(1-exp(-beta*Tinf))*exp(-beta*(time-Tinf))+
                                 (C/gamma)*(1-exp(-gamma*Tinf))*exp(-gamma*(time-Tinf)))) 
   }

   I.Comp3.IVinf.L.optim<-function(x)
   {  logV<-x[1]; 
      logKe<-x[2]; 
      logK12<-x[3];
      logK21<-x[4];
      logK13<-x[5];
      logK31<-x[6];
      ID.list<-unique(tot.data$subject)
      I.value<-0
      for(i in ID.list)
      {  id.s<-which(tot.data$subject==i)
         temp1<- tot.data$conc[id.s]-I.Comp3.IVinf.L(tot.data$dose[id.s],tot.data$Tinf[id.s],
                                       tot.data$time[id.s],logV,logKe,logK12,logK21,logK13,logK31)  
        I.value<-I.value+sum(temp1*temp1)
      }
      I.value
   }

   Three.IVinf.linear<-function(f,f.o,a,b,n.theta,model)
   {  dataMake(a,b)
      setModel(n.theta,model)
      setInitial(f,f.o,model)
   }

   Three.IVinf.linear1<-function(f,f.i,param.list)
   {  Simulation(f,m,K)
      NDfunction<-"numericDeriv(quote(loglike.ind(id.log,T.data,f,logV,logKe,logK12,logK21,logK13,logK31)),c('logV','logKe','logK12','logK21','logK13','logK31'), myenv)"
      FisherInfo(f.i,param.list,NDfunction)
      Pred(f)
      return(tot.data)
   }

  
   ### three Comp zero order absorption   
    Comp3.ZO.L<-function(dose,Tinf,time,mu)
   {  logV<-mu[1]  ; V<-exp(logV)
      logKe<-mu[2] ; Ke<-exp(logKe)
      logK12<-mu[3] ; K12<-exp(logK12)
      logK21<-mu[4] ; K21<-exp(logK21)
      logK13<-mu[5] ; K13<-exp(logK13)
      logK31<-mu[6] ; K31<-exp(logK31)

      a0<-Ke*K21*K31;a1<-Ke*K31+K21*K31+K21*K13+Ke*K21+K31*K12;a2<-Ke+K12+K13+K21+K31
      p<-a1-(a2^2)/3;q<-((2*a2^3)/27)-(a1*a2/3)+a0
      r1<-sqrt(-(p^3/27));r2<-2*r1^(1/3);Phi<-acos(-q/(2*r1))/3
      alpha<--(cos(Phi)*r2-(a2/3))
      beta<--(cos(Phi+(2*pi/3))*r2-(a2/3))
      gamma<--(cos(Phi+(4*pi/3))*r2-(a2/3))

      A<-((K21-alpha)*(K31-alpha))/(V*(alpha-beta)*(alpha-gamma))
      B<-((K21-beta)*(K31-beta))/(V*(beta-alpha)*(beta-gamma))
      C<-((K21-gamma)*(K31-gamma))/(V*(gamma-beta)*(gamma-alpha))
      ifelse(time<=Tk0,
         (dose/Tk0)*(A/alpha*(1-exp(-alpha*time))+B/beta*(1-exp(-beta*time))+C/gamma*(1-exp(-gamma*time))),
         (dose/Tk0)*(A/alpha*(1-exp(-alpha*Tk0))*exp(-alpha*(time-Tk0))+
                      B/beta*(1-exp(-beta*Tk0))*exp(-beta*(time-Tk0))+
                      C/gamma*(1-exp(-gamma*Tk0))*exp(-gamma*(time-Tk0)))) 
   }

   I.Comp3.ZO.L<-function(dose,Tinf,time,logV,logKe,logK12,logK21,logK13,logK31)
   {  V<-exp(logV)
      Ke<-exp(logKe)
      K12<-exp(logK12)
      K21<-exp(logK21)
      K13<-exp(logK13)
      K31<-exp(logK31)

      a0<-Ke*K21*K31;a1<-Ke*K31+K21*K31+K21*K13+Ke*K21+K31*K12;a2<-Ke+K12+K13+K21+K31
      p<-a1-(a2^2)/3;q<-((2*a2^3)/27)-(a1*a2/3)+a0
      r1<-sqrt(-(p^3/27));r2<-2*r1^(1/3);PI<-acos(-q/(2*r1))/3
      alpha<--(cos(PI)*r2-(a2/3))
      beta<--(cos(PI+(2*pi/3))*r2-(a2/3))
      gamma<--(cos(PI+(4*pi/3))*r2-(a2/3))
      A<-((K21-alpha)*(K31-alpha))/(V*(alpha-beta)*(alpha-gamma))
      B<-((K21-beta)*(K31-beta))/(V*(beta-alpha)*(beta-gamma))
      C<-((K21-gamma)*(K31-gamma))/(V*(gamma-beta)*(gamma-alpha))

      A<-((K21-alpha)*(K31-alpha))/(V*(alpha-beta)*(alpha-gamma))
      B<-((K21-beta)*(K31-beta))/(V*(beta-alpha)*(beta-gamma))
      C<-((K21-gamma)*(K31-gamma))/(V*(gamma-beta)*(gamma-alpha))
      ifelse(time<=Tk0,
           (dose/Tk0)*(A/alpha*(1-exp(-alpha*time))+B/beta*(1-exp(-beta*time))+C/gamma*(1-exp(-gamma*time))),
           (dose/Tk0)*(A/alpha*(1-exp(-alpha*Tk0))*exp(-alpha*(time-Tk0))+
                      B/beta*(1-exp(-beta*Tk0))*exp(-beta*(time-Tk0))+
                      C/gamma*(1-exp(-gamma*Tk0))*exp(-gamma*(time-Tk0)))) 
   }

   I.Comp3.ZO.L.optim<-function(x)
   {  logV<-x[1]
      logKe<-x[2]
      logK12<-x[3]
      logK21<-x[4]
      logK13<-x[5]
      logK31<-x[6]
      ID.list<-unique(tot.data$subject)
      I.value<-0
      for(i in ID.list)
      {  id.s<-which(tot.data$subject==i)
         temp1<- tot.data$conc[id.s]-I.Comp3.ZO.L(tot.data$dose[id.s],tot.data$Tinf[id.s],
                                tot.data$time[id.s],logV,logKe,logK12,logK21,logK13,logK31)  
         I.value<-I.value+sum(temp1*temp1)
      }
      I.value
   }

   Three.ZeroOrder.linear<-function(f,f.o,a,b,n.theta,model)
   {  dataMake(a,b)
      setModel(n.theta,model)
      setInitial(f,f.o,model)
   }

   Three.ZeroOrder.linear1<-function(f,f.i,param.list)
   {  Simulation(f,m,K)
      NDfunction<-"numericDeriv(quote(loglike.ind(id.log,T.data,f,logV,logKe,logK12,logK21,logK13,logK31)),
                                        c('logV','logKe','logK12','logK21','logK13','logK31'), myenv)"
      FisherInfo(f.i,param.list,NDfunction)
      Pred(f)
      return(tot.data)
   }
 
   ### plot
 
#### -- TIMEvsDVandPRED.plot
   ID.plot<-function(h,...)
   {  select.data<-Data.Info
      ID<-paste(as.character(sort(unique(select.data$subject))),"   ")
      ID<-matrix(ID,nrow=length(ID))
      colnames(ID)<-c("ID")

      updatePlot<-function(h,...) 
      {  select.data<-Data.Info
         id<-svalue(IDlist)
         data<-select.data[which(select.data$subject==as.numeric(id[1])),]
         par(mfrow=c(1,1))
         plot(data$time,data$conc,type='l',xlab="TIME",ylab="Conc",
            xlim=range(data$time,na.rm=T),
            ylim=range(c(data$conc,data$y.pop,data$y.ind),na.rm=T),
            main=paste("ID:",id[1]))
       
      }
   
      Plot1<-function(h,...) 
      {  par(mfrow=c(1,1))
         select.data<-Data.Info
         plot(select.data$time,select.data$conc,type='n',xlab="TIME",ylab="DV")
         id<-names(table(select.data$subject))
         for(i in 1:length(id))
         {  data<-select.data[which(select.data$subject==id[i]),]
            lines(data$time,data$conc,lty=1,col=1)
          }
      }

      Plot2<-function(h,...) 
      {  select.data<-Data.Info
         ID<-names(table(select.data$subject))
         if(length(ID)>16)       
        {  tkmessageBox(title = "Warning",
                    message = "You have too many individuals. Try [All in one plot]", icon = "error", type = "ok")
        } else
        {
         p<-ceiling(sqrt(length(ID)))
         if(p<=5)
         {  par(mfrow=c(p,p))
            x.r<-c(min(select.data$time,na.rm=T),max(select.data$time,na.rm=T))
            y.r<-c(min(select.data$conc,na.rm=T),max(select.data$conc,na.rm=T))
            for(i in 1:length(ID))
            {   data<-select.data[which(select.data$subject==ID[i]),]
                plot(data$time,data$conc,type='l',xlab="TIME",ylab="Conc",xlim=x.r,ylim=y.r,
                   main=paste("ID:",i))
            }
         }
         }
      }   
    
      IDlist <<- gdroplist(ID,handler=updatePlot) 
      Button1<-gbutton("OK",handler=Plot1)
      Button2<-gbutton("OK",handler=Plot2)
   
      window <- gwindow("ID : Conc vs Time")
      BigGroup <- ggroup(cont=window,anchor=c(-1,1))
      group<-ggroup(cont=BigGroup,horizontal=FALSE)
      tmp<-gframe("ID",container=group)
      add(tmp,IDlist)
      tmp<-gframe("All in one",container=group)
      add(tmp,Button1)
      tmp<-gframe("Individual groups",container=group)
      add(tmp,Button2,expand=TRUE)
      tmp<-gframe("DV   : black",container=group)
      tmp<-gframe("pop : blue",container=group)
      tmp<-gframe("IND : red",container=group)  
      add(BigGroup, ggraphics())
   }
#### -- DVvsPRED.plot
   DVvsPRED.plot<-function(h,...)
   {  updatePlot<-function(h,...)
      {  condX.V <-svalue(VarList.X)
         condY.V<-svalue(VarList.Y)
         select.data<-OUTPUT.tot[[RUN.id]]$y.pred
         if(!is.na(condX.V)& !is.na(condY.V))
         {  X<-select.data$conc
            Y<-select.data[,condY.V]
            plot(X,Y,xlab=condX.V,ylab=condY.V)
            abline(a=0,b=1,col=2)
         }
      }
      VarList.X<-gdroplist("conc")
      VarList.Y<-gdroplist(c("y.pop","y.ind"))
      Button<-gbutton("OK",handler=updatePlot)
    
      win<-gwindow("DV vs PRED plot")
      BigGroup<-ggroup(cont=win)
      group=ggroup(horizontal=FALSE,cont=BigGroup)
      tmp=gframe(" X variable",container=group)
      add(tmp,VarList.X)
      tmp=gframe(" Y variable",container=group)
      add(tmp,VarList.Y)

      tmp=gframe("Plot",container=group)
      add(tmp,Button,expand=TRUE)
      add(BigGroup,ggraphics())
   }

#### -- DVvsRES.plot
   DVvsRES.plot<-function(h,...)
   {  updatePlot<-function(h,...)
      {  par(mfrow=c(1,1))
         condX.V <-svalue(VarList.X)
         condY.V<-svalue(VarList.Y)
         select.data<-OUTPUT.tot[[RUN.id]]$y.pred
         if(!is.na(condX.V)& !is.na(condY.V))
         {  X<-select.data$conc
            Y<-select.data[,condY.V]
            plot(X,Y,xlab=condX.V,ylab=condY.V,ylim=c(-4,4))
            abline(h=0,col=3,lwd=1.5)
            if(condY.V=="y.ind.wres")
            {  abline(h=2,col=2,lwd=1.5)
               abline(h=-2,col=2,lwd=1.5)
            }

         }
      }

      Plot1<-function(h,...) 
      {  par(mfrow=c(1,2))
         plot(OUTPUT.tot[[RUN.id]]$y.pred[,"conc"],
              OUTPUT.tot[[RUN.id]]$y.pred[,"y.ind.res"],xlab="DV",ylab="RES",ylim=c(-4,4))
            abline(h=0,col=3,lwd=1.5)
         plot(OUTPUT.tot[[RUN.id]]$y.pred[,"conc"],
             OUTPUT.tot[[RUN.id]]$y.pred[,"y.ind.wres"],xlab="DV",ylab="WRES",ylim=c(-4,4))
            abline(h=0,col=3,lwd=1.5)
            abline(h=2,col=2,lwd=1.5)
            abline(h=-2,col=2,lwd=1.5)
      }
    
      VarList.Y<-gdroplist(c("y.ind.res","y.ind.wres"))
  
      VarList.X<-gdroplist(c("conc"))  
      Button1<-gbutton("OK",handler=updatePlot)
      Button2<-gbutton("OK",handler=Plot1)
    
      win<-gwindow("DV vs RES plot")
      BigGroup<-ggroup(cont=win)
      group<-ggroup(horizontal=FALSE,cont=BigGroup)
      tmp<-gframe(" X variable",container=group)
      add(tmp,VarList.X)
      tmp<-gframe(" Y variable",container=group)
      add(tmp,VarList.Y)
      tmp<-gframe("Plot",container=group)
      add(tmp,Button1,expand=TRUE)
      tmp<-gframe("All Plots",container=group)
      add(tmp,Button2,expand=TRUE) 
      add(BigGroup,ggraphics())
   }

#### -- TIMEvsRES.plot
   TIMEvsRES.plot<-function(h,...)
   {  updatePlot<-function(h,...)
      {  condX.V <-svalue(VarList.X)
         condY.V<-svalue(VarList.Y)
         select.data<-OUTPUT.tot[[RUN.id]]$y.pred
         if(!is.na(condX.V)& !is.na(condY.V))
         {  X<-select.data$time
            Y<-select.data[,condY.V]
            plot(X,Y,xlab=condX.V,ylab=condY.V,type='n')
            id<-as.numeric(names(table(select.data$subject)))
            for(i in 1:length(id))
            {  X.data<-X[which(select.data$subject==id[i])]
               Y.data<-Y[which(select.data$subject==id[i])]
               lines(X.data,Y.data,lty=2)
            }
            abline(h=0,col=2)
         }
      }

      VarList.Y<-gdroplist(c("y.ind.res","y.ind.wres"))
   
      VarList.X<-gdroplist(c("time"))
      Button<-gbutton("OK",handler=updatePlot)
    
      win<-gwindow("TIME vs RES plot")
      BigGroup<-ggroup(cont=win)
      group<-ggroup(horizontal=FALSE,cont=BigGroup)
      tmp<-gframe(" X variable",container=group)
      add(tmp,VarList.X)
      tmp<-gframe(" Y variable",container=group)
      add(tmp,VarList.Y)
      tmp<-gframe("Plot",container=group)
      add(tmp,Button,expand=TRUE)
      add(BigGroup,ggraphics())
   }

#### -- TIMEvsDVandPRED.plot
   TIMEvsDVandPRED.plot<-function(h,...)
   {  select.data<-OUTPUT.tot[[RUN.id]]$y.pred
      ID<-paste(as.character(sort(unique(select.data$subject))),"   ")
      ID<-matrix(ID,nrow=length(ID))
      colnames(ID)<-c("ID")

      updatePlot<-function(h,...) 
      {  select.data<-OUTPUT.tot[[RUN.id]]$y.pred
         id<-svalue(IDlist)
         data<-select.data[which(select.data$subject==as.numeric(id[1])),]
         par(mfrow=c(1,1))
         plot(data$time,data$conc,type='l',xlab="TIME",ylab="Conc",
            xlim=range(data$time,na.rm=T),
            ylim=range(c(data$conc,data$y.pop,data$y.ind),na.rm=T),
            main=paste("ID:",id[1]))
         lines(data$time,data$y.pop,lty=1,col=4)
         lines(data$time,data$y.ind,lty=1,col=2)       
      }
   
      Plot1<-function(h,...) 
      {  par(mfrow=c(1,1))
         select.data<-OUTPUT.tot[[RUN.id]]$y.pred
         plot(select.data$time,select.data$conc,type='n',xlab="TIME",ylab="DV")
         id<-names(table(select.data$subject))
         for(i in 1:length(id))
         {  data<-select.data[which(select.data$subject==id[i]),]
            lines(data$time,data$conc,lty=1,col=1)
            lines(data$time,data$y.pop,lty=1,col=4)
            lines(data$time,data$y.ind,lty=1,col=2)
         }
      }

      Plot2<-function(h,...) 
      {  select.data<-OUTPUT.tot[[RUN.id]]$y.pred
         ID<-names(table(select.data$subject))
         if(length(ID)>16)       
        {  tkmessageBox(title = "Warning",
                    message = "You have too many individuals. Try [All in one plot]", icon = "error", type = "ok")
        } else
        {
         p<-ceiling(sqrt(length(ID)))
         if(p<=5)
         {  par(mfrow=c(p,p))
            x.r<-c(min(select.data$time,na.rm=T),max(select.data$time,na.rm=T))
            y.r<-c(min(select.data$conc,na.rm=T),max(select.data$conc,na.rm=T))
            for(i in 1:length(ID))
            {   data<-select.data[which(select.data$subject==ID[i]),]
                plot(data$time,data$conc,type='l',xlab="TIME",ylab="Conc",xlim=x.r,ylim=y.r,
                   main=paste("ID:",i))
                lines(data$time,data$y.pop,lty=1,col=4)
                lines(data$time,data$y.ind,lty=1,col=2)
            }
         }
         }
      }   
    
      IDlist <<- gdroplist(ID,handler=updatePlot) 
      Button1<-gbutton("OK",handler=Plot1)
      Button2<-gbutton("OK",handler=Plot2)
   
      window <- gwindow("ID : Conc vs Time")
      BigGroup <- ggroup(cont=window,anchor=c(-1,1))
      group<-ggroup(cont=BigGroup,horizontal=FALSE)
      tmp<-gframe("ID",container=group)
      add(tmp,IDlist)
      tmp<-gframe("All in one",container=group)
      add(tmp,Button1)
      tmp<-gframe("Individual groups",container=group)
      add(tmp,Button2,expand=TRUE)
      tmp<-gframe("DV   : black",container=group)
      tmp<-gframe("pop : blue",container=group)
      tmp<-gframe("IND : red",container=group)  
      add(BigGroup, ggraphics())
   }

#### -- EBEvsCOV.plot1
   EBEvsCOV.plot<-function(h,...)
   {  ETA.list<-NULL
      N.eta<-OUTPUT.tot[[RUN.id]]$n.theta
      for(i in 1:N.eta)
         ETA.list<-c(ETA.list,paste("ETA(",i,")",sep=""))
      id.list<-unique(Data.Info$subject)
      ETA<-tot.data$theta.est
      if(length( which(Var.Prop=="COV"))!=0)
     {   temp.id<-Data.Info[,which(Var.Prop=="COV")[1]]
          if(length(unique( temp.id[Data.Info[,"subject"]==id.list[1]]))!=1)
          {   COV.list<-c(which(Var.Prop=="ID"),which(Var.Prop=="TIME"),which(Var.Prop=="COV"))
               select.data<-as.matrix(Data.Info[,COV.list])
          } else
          {   COV.list<-c(which(Var.Prop=="ID"),which(Var.Prop=="COV"))
               select.data<-as.matrix(Data.Info[,COV.list])
               select.data<-select.data[!duplicated(select.data[,1]),]
          } 
          COV.list<-(which(Var.Prop=="COV"))
          Var.Name<-colnames(Data.Info)
 
          updatePlot<-function(h,...)
          {  condX.V <-svalue(VarList.X,index=T)
             condY.V<-svalue(VarList.Y,index=T)

             if(!is.na(condX.V)& !is.na(condY.V))
             {  Y<<-select.data[,condY.V+1]
                 X<<-OUTPUT.tot[[RUN.id]]$theta[,condX.V]
                plot(X,Y,xlab=ETA.list[condX.V],ylab=Var.Name[COV.list[condY.V]])
             }
          }
          linearD<-function(h,...)
         { 
             if(length(table(Y))==2)
             {   lines(sort(X),sort(glm(Y~X,family=binomial)$fitted),lwd=2,col=2)
             } else
             {   abline(lm(Y~X),lwd=2,col=2)
             } 
          }
           lowessD<-function(h,...)
         {  lines(lowess(Y~X),lwd=2,col=2)
          }

          VarList.Y<-gdroplist(Var.Name[COV.list])
          VarList.X<-gdroplist(ETA.list)

          Button1<-gbutton("OK",handler=updatePlot)
          Button2<-gbutton("linear",handler=linearD)
          Button3<-gbutton("lowess",handler=lowessD)

    
          win<-gwindow("XY plot")
          BigGroup<-ggroup(cont=win)
          group<-ggroup(horizontal=FALSE,cont=BigGroup)
          tmp<-gframe(" X variable",container=group)
          add(tmp,VarList.X)
          tmp<-gframe(" Y variable",container=group)
          add(tmp,VarList.Y)

          tmp<-gframe("Plot",container=group)
          add(tmp,Button1,expand=TRUE)
          tmp<-gframe("lines",container=group)
          add(tmp,Button2,expand=TRUE)
          add(tmp,Button3,expand=TRUE)
          add(BigGroup,ggraphics())
    }
}

  
   ### simulation
   inv.m<-function(A)
   {  temp.eig<-eigen(A)
      temp.eig$vectors%*%diag(1/temp.eig$values)%*%t(temp.eig$vectors)
   }

   Phipost.ind<-function(f,theta,data,data.t)
   {  a<-data.t$a; b<-data.t$b;c<-data.t$c
      f.value<-f(data$dose,data$Tinf,data$time,theta)
      f.value[f.value<0.0000001]<-0
      A1<- -log(a+b*f.value^c)
      A2<--0.5*(((data$conc-f.value)/(a+b*f.value^c))^2)
      A3<--0.5*matrix((theta-data.t$mu),ncol=length(theta))%*%solve(data.t$omega)%*%
                                     matrix((theta-data.t$mu),nrow=length(theta))
      return(sum(A1)+sum(A2)+sum(A3))
   }

   logmultinorm <- function(x, m, v) 
   {  return(-0.5 * t(x - m) %*% solve(v) %*% (x - m))
   }

   dataMake<-function(a,b)
   {  tot.data<-list()
      tot.data$Dose.info<-Dose.Info
      tot.data$Data.info<-Data.Info
      tot.data$conc<-Data.Info$conc
      tot.data$time<-Data.Info$time
      tot.data$subject<-Data.Info$subject
      id.list<-unique(tot.data$subject)
      tot.data$dose<-rep(NA,length(tot.data$conc))
      tot.data$Tinf<-rep(NA,length(tot.data$conc))

      for(i in id.list)
      {  tot.data$dose[tot.data$subject==i]<-Dose.Info$dose[Dose.Info$subject==i]
         if(!is.null(Dose.Info$Tinf))
            tot.data$Tinf[tot.data$subject==i]<-Dose.Info$Tinf[Dose.Info$subject==i]
      }
   
      tot.data$a.f<-T
      tot.data$b.f<-T
      tot.data$c.f<-T

      tot.data$N<-length(table(tot.data$subject))
      tot.data<<-tot.data
   }

   setModel<-function(n.theta,model)
   {  temp<-tot.data
      temp$n.theta<-n.theta
      temp$model<-model
      tot.data<<-temp
   }

   setInitial<-function(f,f.o,PKmethod)
   {  id.list<-tot.data$id.list
      N<-tot.data$N
      temp.tot.data<-tot.data
      if( admin=="First Order absorption")
      {  init.value<-c(3,0,-2)
         if(tot.data$n.theta>3)
         {  init.value<-c(init.value,rep(0.1,tot.data$n.theta-3))
         }
      } else
      {  init.value<-c(3,-2)
         if(tot.data$n.theta>2)
         {  init.value<-c(init.value,rep(0,tot.data$n.theta-2))
         }
      } 
      mu.temp<-optim(init.value,f.o,method="L-BFGS-B")$par
      
      initial.GUI<-gwindow("Initial Values",width=50)
      g<-ggroup(horizontal=FALSE,cont=initial.GUI)
      a<-list()
      for(i in 1:length(param.list))
           a[[i]]<-gedit("0")
      for(j in 1:length(param.list))
          svalue(a[[j]])<-round(mu.temp[j],2)
      button<-gbutton("OK",
              handler=function(h,...)
              {  temp.tot.data<-tot.data
                 mu.keep<-NULL
                 for(j in 1:length(param.list))
                    mu.keep<-c(mu.keep,svalue(a[[j]]))
                 dispose(initial.GUI)   
                 mu.keep<<-as.numeric(mu.keep)
                 temp.tot.data$mu<-as.numeric(mu.keep)
                 temp.tot.data$theta<-NULL
                 for(i in 1:tot.data$N)
                    temp.tot.data$theta<-rbind(temp.tot.data$theta,c(temp.tot.data$mu))
                 temp.tot.data$omega<-diag(rep(0.1,each=temp.tot.data$n.theta))
                 temp.tot.data$a<-0.1
                 temp.tot.data$b<-0.1
                 temp.tot.data$c<-0.1
                 s1<-temp.tot.data$theta
                 s2<-t(temp.tot.data$theta)%*%temp.tot.data$theta
                 diff.sq<-NULL
                 for(id.i in id.list)
                 {   temp.list<-which(temp.tot.data$subject==id.i)
                     f.value<-f(temp.tot.data$dose[temp.list],temp.tot.data$Tinf[temp.list],
                                temp.tot.data$time[temp.list], temp.tot.data$theta[which(id.list==id.i),])
                     temp.sq<-(f.value-temp.tot.data$conc[temp.list])^2
                     if(temp.tot.data$b.f)
                     {   temp.sq[f.value!=0]<-temp.sq[f.value!=0]/f.value[f.value!=0]^2
                         temp.sq[f.value==0]<-0
                     }                                                         
                     diff.sq<-c(diff.sq,temp.sq)
                 }
                 s3<-sum(diff.sq)
                 temp<-list()
                 temp$s1<-s1; temp$s2<-s2;temp$s3<-s3
                 Init.v<<-temp
                 tot.data<<-temp.tot.data
                 if(comp=="one"&admin=="First Order absorption"&elim=="linear")
                 {   result<-One.FirstOrder.linear1(Comp1.FO.L,I.Comp1.FO.L,param.list)
                 } else if(comp=="one"&admin=="IV bolus"&elim=="linear")
                 {   result<-One.IVbolus.linear1(Comp1.IVbolus.L,I.Comp1.IVbolus.L,param.list)
                 } else if(comp=="one"&admin=="IV infusion"&elim=="linear")
                 {   result<-One.IVinf.linear1(Comp1.IVinf.L,I.Comp1.IVinf.L,param.list)
                 } else if(comp=="one"&admin=="Zero Order absorption"&elim=="linear")
                 {   result<-One.ZeroOrder.linear1(Comp1.ZO.L,I.Comp1.ZO.L,param.list)
                 } else if(comp=="two"&admin=="First Order absorption"&elim=="linear")
                 {   result<-Two.FirstOrder.linear1(Comp2.FO.L,I.Comp2.FO.L,param.list)
                 } else if(comp=="two"&admin=="IV bolus"&elim=="linear")
                 {   result<-Two.IVbolus.linear1(Comp2.IVbolus.L,I.Comp2.IVbolus.L,param.list)
                 } else if(comp=="two"&admin=="IV infusion"&elim=="linear")
                 {   result<-Two.IVinf.linear1(Comp2.IVinf.L,I.Comp2.IVinf.L,param.list)
                 } else if(comp=="two"&admin=="Zero Order absorption"&elim=="linear")
                 {   result<-Two.ZeroOrder.linear1(Comp2.ZO.L,I.Comp2.ZO.L,param.list)
                 } else if(comp=="three"&admin=="First Order absorption"&elim=="linear")
                 {   result<-Three.FirstOrder.linear1(Comp3.FO.L,I.Comp3.FO.L,param.list)
                 } else if(comp=="three"&admin=="IV bolus"&elim=="linear")
                 {   result<-Three.IVbolus.linear1(Comp3.IVbolus.L,I.Comp3.IVbolus.L,param.list)
                 } else if(comp=="three"&admin=="IV infusion"&elim=="linear")
                 {   result<-Three.IVinf.linear1(Comp3.IVinf.L,I.Comp3.IVinf.L,param.list)
                 } else if(comp=="three"&admin=="Zero Order absorption"&elim=="linear")
                 {   result<-Three.ZeroOrder.linear1(Comp3.ZO.L,I.Comp3.ZO.L,param.list)
                 }  
                 TextMake(model,param.list,result)
  	             temp<-OUTPUT.tot 
     	           temp[[RUN.id]]<-result
     	           OUTPUT.tot<<-temp
     	           output.window<-gwindow(paste("Estimation Result :", RUN.id))
     	           g<-ggroup(horizontal=FALSE,cont=output.window)
     	           tmp<-gframe("Output",container=g)
     	           a<-gtext(output.txt1,width=650,height=500,
                                 font.attr=c(sizes="large",family="monospace"))
     	          add(tmp,a)
              })
      for(i in 1:length(param.list))
      {  tmp<-gframe(param.list[i],cont=g)
         add(tmp,a[[i]])
      }
      tmp<-gframe("",cont=g)
      add(tmp,button)
   }

   funMakeN<-function()
   {  LST<-paste(      "loglike.ind.N <-function(id.log,T.data,f) " ,sep="") 
      LST<-paste(LST,   "{   id.list<-unique(T.data$subject) ;")
      LST<-paste(LST,   "    sel.id<-which(T.data$subject==id.log) ;")
      LST<-paste(LST,          " conc<-T.data$conc[sel.id]; ")
      LST<-paste(LST,          " Tinf<-T.data$Tinf[sel.id];")
      LST<-paste(LST,         "  time<-T.data$time[sel.id];")
      LST<-paste(LST,         "  dose<-T.data$dose[sel.id];")
      LST<-paste(LST,        "  ni<-length(conc);")
      LST<-paste(LST,         "  a<-T.data$a;")
      LST<-paste(LST,         "  b<-T.data$b;")
      LST<-paste(LST,         "  c<-T.data$c;")
      LST<-paste(LST,        "  mu<-matrix(T.data$mu,nrow=T.data$n.theta);")
      LST<-paste(LST,         "  abf<-a+b*f(dose,Tinf,time,mu)^c ;")
      LST<-paste(LST,         "  abf<-ifelse(abf==0,1,abf);")
      LST<-paste(LST,         "  omega<-T.data$omega;")
      LST<-paste(LST,         "  theta.i<-T.data$theta[which(id.list==id.log),];")
      LST<-paste(LST,         "  -0.5*(1+ni)*log(2*pi) -sum(log(abf^2))- ")
      LST<-paste(LST,         "       0.5*sum((conc-  f(dose,Tinf,time,mu))^2/(abf)^2)-")
      LST<-paste(LST,         "       0.5*log(det(omega))-0.5*t((theta.i-mu))%*%solve(omega)%*%(theta.i-mu);}")
      write(LST,"temp.txt")
      source("temp.txt")
   }

   Simulation<-function(f,m,K)
   {  id.list<-unique(tot.data$subject)
      temp.tot.data<-tot.data
      temp.v<-Init.v
      temp.data<-list()
      temp.data$mu<-temp.tot.data$mu
      temp.data$omega<-temp.tot.data$omega 
      temp.data$sigma<-temp.tot.data$sigma
      temp.data$theta<-temp.tot.data$theta
      temp.data$N<-tot.data$N
      theta.TOT<-list()
      KEEP.tot<-list()
      theta.i<-matrix(0,ncol=(tot.data$n.theta+1),nrow=tot.data$N)
      omega.i<-matrix(0,nrow=tot.data$n.theta,ncol=tot.data$n.theta)
      mu.i.2<-matrix(0,nrow=tot.data$n.theta)
      funMakeN()
      L.keep<-NULL
      temp.tot.data$L<--10^10

      for(k in 1:K)
      {  temp.old<-temp.tot.data
         s1.old<-temp.v$s1;s2.old<-temp.v$s2;s3.old<-temp.v$s3
         temp.data$subject<-NULL
         theta.start<-temp.data$mu
         theta.all<-NULL
         theta.se<-NULL
         f.se<-NULL
         i<-0
         for(id.i in id.list)
         {  i<-i+1
            temp.data$conc<-temp.tot.data$conc[temp.tot.data$subject==id.i]
            temp.data$time<-temp.tot.data$time[temp.tot.data$subject==id.i]
            temp.data$dose<-ifelse(length(temp.tot.data$dose)==0,NA,
                                temp.tot.data$dose[temp.tot.data$subject==id.i])
            temp.data$Tinf<-ifelse(length(temp.tot.data$Tinf)==0,NA,
                                 temp.tot.data$Tinf[temp.tot.data$subject==id.i])
            th0<-theta.start
            pb <- length(th0)
            mu<-theta.start
            cv <-  chol(temp.tot.data$omega)
            f0 <- Phipost.ind(f,th0, temp.data,temp.tot.data)
            KEEP<-NULL
            f.value<-NULL
            for(j in 1:m) 
            {  th1 <- th0 + t(cv) %*% array(rnorm(pb), c(pb, 1))
               f1 <- Phipost.ind(f,th1,temp.data,temp.tot.data)
               R <- exp( f1 - f0)
               u <- runif(1) < R
               if(u == 1) 
               {  th0 <- th1
                  f0 <- f1
                  KEEP<-rbind(KEEP,c(th0))
                  temp1<-f(temp.data$dose,temp.data$Tinf,temp.data$time,  th0)
                  temp2<-temp.data$conc-temp1
                  f.value<-rbind(f.value,cbind(temp1,temp2))
               }
            }
            f.group<-rep(1:length(temp.data$time),nrow(KEEP))
            f.se<-c(f.se,tapply(f.value[,2],f.group,sd))
            KEEP.n<-min(nrow(KEEP),10)
            KEEP<-KEEP[(nrow(KEEP)-KEEP.n+1) : nrow(KEEP),]
            if(KEEP.n==1)
               KEEP<-rbind(KEEP,KEEP)
            KEEP.tot[[id.i]]<-KEEP
            f.value<-f.value[(nrow(f.value)-length(temp.data$time)*KEEP.n+1) : nrow(f.value),]
            theta.i[i,]<-c(id.i,apply(KEEP,2,mean))
            theta.se<-rbind(theta.se,apply(KEEP,2,sd))
            mu.i.2<-mu.i.2+theta.i[i,-1,drop=TRUE]
         }
         mu<-apply(theta.i[,-1],2,mean)
         omega.i<-matrix(0,nrow=tot.data$n.theta,ncol=tot.data$n.theta) 
         for(i in 1:tot.data$N)
            omega.i<-omega.i+(theta.i[i,-1]-mu)%*%t(theta.i[i,-1]-mu)
         omega<-omega.i/tot.data$N
         fr<-function(x)
         {  a<-x[1]
            b<-x[2]
            c<-x[3]
            fij<- f(tot.data$dose,tot.data$Tinf,tot.data$time,temp.tot.data$mu)
            fij[fij<0.000001]<-0
            sum(log(a+b*fij^c)+0.5*((tot.data$conc-fij)/(a+b*fij^c))^2)
         }
         abc<-optim(c(temp.tot.data$a,temp.tot.data$b,temp.tot.data$c),fr)$par
         a<-abc[1]
         b<-abc[2]
         c<-abc[3]
         L<-0
 
         temp.tot.data$mu<-mu
         temp.tot.data$omega<-omega
         temp.tot.data$a<-a
         temp.tot.data$b<-b
         temp.tot.data$c<-c
         temp.tot.data$theta<-theta.i[,-1]
         temp.tot.data$theta.se<-theta.se

         for(id.i in id.list)
         {  L<- L+loglike.ind.N(id.i,temp.tot.data,f)  
         }
         temp.tot.data$L<-L
         L.keep<-rbind(L.keep,c(k,L))
         tot.data<<-temp.tot.data
         L.KEEP<<-L.keep  
      }
      tot.data<<-temp.tot.data
      L.KEEP<<-L.keep
   }

   funMake<-function(paramlist)
   {  ftn.name<-"f.i(dose,Tinf,time"
      param.name<-NULL
      for(i in paramlist)
      {  if(i==paramlist[1])
         {   param.name<- paste(param.name,as.character(i))
         } else
         {   param.name<- paste(param.name,as.character(i),sep=",")
         }
         ftn.name<-paste(ftn.name,as.character(i),sep=",")
      }
      ftn.name<-paste(ftn.name,")")

      LST<-paste(      "loglike.ind <-function(id.log,T.data,f,",param.name,") " ,sep="") 
      LST<-paste(LST,   "{   id.list<-unique(T.data$subject) ;")
      LST<-paste(LST,   "    sel.id<-which(T.data$subject==id.log) ;")
      LST<-paste(LST,          " conc<-T.data$conc[sel.id]; ")
      LST<-paste(LST,          " Tinf<-T.data$Tinf[sel.id];")
      LST<-paste(LST,         "  time<-T.data$time[sel.id];")
      LST<-paste(LST,         "  dose<-T.data$dose[sel.id];")
      LST<-paste(LST,        "  ni<-length(conc);")
      LST<-paste(LST,         "  a<-T.data$a;")
      LST<-paste(LST,         "  b<-T.data$b;")
      LST<-paste(LST,         "  c<-T.data$c;")
      LST<-paste(LST,         "  abf1<-f(dose,Tinf,time,", param.name,");")
      LST<-paste(LST,         "  abf<-a+b*ifelse(abf1<0.000001,0,abf1)^c ;")
      LST<-paste(LST,         "  abf<-ifelse(abf==0,1,abf);")
      LST<-paste(LST,        "  mu<-matrix(T.data$mu,nrow=T.data$n.theta);")
      LST<-paste(LST,         "  omega<-T.data$omega;")
      LST<-paste(LST,         "  theta.i<-T.data$theta[which(id.list==id.log),];")
      LST<-paste(LST,         "  -0.5*(1+ni)*log(2*pi) - sum(log(abf))-")
      LST<-paste(LST,         "       0.5*sum((conc-  f(dose,Tinf,time,",param.name,"))^2/(abf)^2)-")
      LST<-paste(LST,         "       0.5*log(det(omega))-0.5*t((theta.i-mu))%*%solve(omega)%*%(theta.i-mu);}")
      write(LST,"temp.txt")
      source("temp.txt")
   }

   FisherInfo<-function(f.i,param.list,NDfunction)
   {  id.list<-unique(tot.data$subject)
      temp.tot.data<-tot.data
      mu<-tot.data$mu
      N<-tot.data$N
      myenv<-NULL
      myenv<-new.env() 
      for(i in 1:length(param.list))
      {  assign(param.list[i],mu[i],envir=myenv)
      }
      assign("T.data", tot.data, envir = myenv)
      assign("f",f.i , envir = myenv)
      Keep.env<<-myenv
      Jacob.save<-array(0, c(N, length(mu)))

      for(id.i in id.list)
      {  funMake(param.list)
         assign("id.log",id.i , envir = myenv)
         temp<-eval(parse(text=NDfunction))
         Jacob.save[which(id.list==id.i),]<-attr(temp,"gradient")
      }
      temp<-t(Jacob.save)%*%Jacob.save

      temp.tot.data$mu.var<-inv.m(temp)
      tot.data<<-temp.tot.data
   }

   Pred<-function(f)
   {  temp.tot.data<-tot.data
      id.list<-unique(temp.tot.data$subject)
      temp.data<-list()
      temp.data$mu<-temp.tot.data$mu
      temp.data$omega<-temp.tot.data$omega 
      temp.data$sigma<-temp.tot.data$sigma
      temp.data$theta<-temp.tot.data$theta

      y.pred<-NULL
      for(id.i in id.list)
      {  temp.data$dose<-tot.data$dose[tot.data$subject==id.i]
         temp.data$Tinf<-tot.data$Tinf[tot.data$subject==id.i]
         temp.data$conc<-tot.data$conc[tot.data$subject==id.i]
         temp.data$time<-tot.data$time[tot.data$subject==id.i]
         theta.i<-temp.data$theta[which(id.list==id.i),]
         temp<-cbind(id.i,temp.data$time,temp.data$conc,f(temp.data$dose,temp.data$Tinf,
                  temp.data$time,temp.data$mu),
         f(temp.data$dose,temp.data$Tinf,temp.data$time,theta.i))
         y.pred<-rbind(y.pred,temp)
      }  
      colnames(y.pred)<-c("subject","time","conc","y.pop","y.ind")
      y.pred<-data.frame(y.pred)
      y.pop.res<-y.pred$conc-y.pred$y.pop
      y.pop.wres<-y.pop.res/temp.tot.data$f.se
      y.ind.res<-y.pred$conc-y.pred$y.ind
      y.ind.res[y.pred$conc==0]<-0
      s.ind <-(tot.data$a+tot.data$b*y.pred$y.ind^tot.data$c)^2
      y.ind.wres<-y.ind.res/sqrt(s.ind)
      y.pred$y.pop.res<-y.pop.res
      y.pred$y.pop.wres<-y.pop.wres
      y.pred$y.ind.res<-y.ind.res
      y.pred$y.ind.wres<-y.ind.wres
      temp.tot.data$y.pred<-y.pred
      temp.tot.data$MSE<-sum((tot.data$conc-y.pred$y.ind)^2)/length(tot.data$conc)
      tot.data<<-temp.tot.data
   }
   ### GUI related
   MakeDose<-function(D.temp)
   {  Var.Name<-colnames(D.temp)
      var.list<-c("ID","time","dose","amt","rate")
      a<-list()
      for(i in 1:length(Var.Name))
          a[[i]]<-gdroplist(var.list)
      button<-gbutton("OK",handler=function(h,...)
              {  Var.Prop.Keep<-NULL
                 for(i in 1:length(Var.Name))
                    Var.Prop.Keep<-c(Var.Prop.Keep,svalue(a[[i]]))
                 colnames(D.temp)[which(Var.Prop.Keep=="ID")]<-"subject"
                 colnames(D.temp)[which(Var.Prop.Keep=="time")]<-"time"
                 colnames(D.temp)[which(Var.Prop.Keep=="dose")]<-"dose"
                 colnames(D.temp)[which(Var.Prop.Keep=="amt")]<-"amt"
                 colnames(D.temp)[which(Var.Prop.Keep=="rate")]<-"rate"
                 if(!is.null(D.temp$amt)& !is.null(D.temp$rate))
                 {  D.temp$Tinf<-D.temp$amt/D.temp$rate
                    D.temp$dose<-D.temp$amt
                 }
                 Dose.Info<<-D.temp
                 dispose(tt)
               })
      tt <- gwindow("Choose Variable Properties",width=50)
      g<-ggroup(horizontal=FALSE,cont=tt)
      tmp<-gframe("Choose Variable Properties",cont=g)
      for(i in 1:length(Var.Name))
      {  tmp<-gframe(Var.Name[i],cont=g)   
         add(tmp,a[[i]])
      }
      tmp<-gframe("",cont=g)   
      add(tmp,button)
   }


   MakeDV<-function(D.temp)
   {  Var.Name<-colnames(D.temp)
      var.list<-c("COV","ID","TIME","DV","Ignore                  ")
      a<-list()
      for(i in 1:length(Var.Name))
         a[[i]]<-gdroplist(var.list)
      button<-gbutton("OK",handler=function(h,...)
              {   Var.Prop.Keep<-NULL
                  for(i in 1:length(Var.Name))
                     Var.Prop.Keep<-c(Var.Prop.Keep,svalue(a[[i]]))
                  colnames(D.temp)[which(Var.Prop.Keep=="ID")]<-"subject"
                  colnames(D.temp)[which(Var.Prop.Keep=="DV")]<-"conc"
                  colnames(D.temp)[which(Var.Prop.Keep=="TIME")]<-"time"
                  Data.Info<<-D.temp
                  Var.Prop<<-Var.Prop.Keep
                  dispose(tt)
              })
      tt <- gwindow("Choose Variable Properties",width=50)
      g<-ggroup(horizontal=FALSE,cont=tt)
      tmp<-gframe("Choose Variable Properties",cont=g)
      for(i in 1:length(Var.Name))
      {  tmp<-gframe(Var.Name[i],cont=g)   
         add(tmp,a[[i]])
      }
      tmp<-gframe("",cont=g)   
      add(tmp,button)
   }


   TextMake<-function(model,param.list,result)
   {   mu.result<-round(result$mu,3)
       mu.result.exp<-round(exp(result$mu),3)
       mu.se.result<-round(sqrt(diag(result$mu.var)),5)
       sigma.result.a<-round(result$a,5)
       sigma.result.b<-round(result$b,5)
       sigma.result.c<-round(result$c,5)
       result.MSE<-round(result$MSE,5)
       omega.result<-round(result$omega,4)
       output.txt1<-paste("--------------------------------------------------------------------------------\n")
       output.txt1<-paste(output.txt1,model,"- combined error\n")
       output.txt1<-paste(output.txt1,"--------------------------------------------------------------------------------\n")
       output.txt1<-paste(output.txt1,"MSE : ",result.MSE,"\n")
       output.txt1<-paste(output.txt1,"--------------------------------------------------------------------------------\n")
       output.txt1<-paste(output.txt1,"THETA\n")
       for( i in 1:length(param.list))
       {  output.txt1<-paste(output.txt1,param.list[i],"  : ",mu.result[i],
               " (", mu.se.result[i],") \t", strsplit(param.list[i],split="log")[[1]][2],":",mu.result.exp[i],"\n")
       }
       output.txt1<-paste(output.txt1,"SIGMA\n")
       output.txt1<-paste(output.txt1,"     ",sigma.result.a,"+", sigma.result.b,"* F^", sigma.result.c,"\n")
       output.txt1<-paste(output.txt1,"OMEGA\n")
       for(i in 1:nrow(omega.result))
       {   for(j in 1:nrow(omega.result)) 	
               output.txt1<-paste(output.txt1,"     ",omega.result[i,j])
           output.txt1<-paste(output.txt1,"\n")
       }
       output.txt1<-paste(output.txt1,"--------------------------------------------------------------------------------\n")
       output.txt1<<-output.txt1
   }

   PrintResult<-function()
   {   TextMake(model,param.list,result)
       temp<-OUTPUT.tot 
       temp[[RUN.id]]<-result
       OUTPUT.tot<<-temp
       output.window<-gwindow(paste("Estimation Result :", RUN.id))
       g<-ggroup(horizontal=FALSE,cont=output.window)
       tmp<-gframe("Output",container=g)
       a<-gtext(output.txt1,width=650,height=500,
               font.attr=c(sizes="large",family="monospace"))
       add(tmp,a)
   }
    
   ### GUI 	
    	
   OpenDoseFile<-function(h,...)
   {   file.dose<<-tclvalue(tkgetOpenFile(title="Open Dosing Information File"))
       D.temp<-read.csv(file.dose,na.strings=".")
       MakeDose(D.temp)
   }

   OpenDataFile<-function(h,...)
   {   file.dose<<-tclvalue(tkgetOpenFile(title="Open Data File"))
       D.temp<-read.csv(file.dose,na.strings=".")
       MakeDV(D.temp)
   }

	 runPK<-function(h,...)
	 { 	RUN.id<<- RUN.id+1
	   	comp<<-svalue(level1)
		  admin<<-svalue(level2)
   		elim<<-svalue(level3)
   		error<<-svalue(level4)
      K<<-svalue( KWidget)
      m<<-svalue( mWidget)
 		  a<-1;b<-1
 
      ## One-comp / First Order Abroption
   		if(comp=="one"&admin=="First Order absorption"&elim=="linear")
   		{	  if(!is.null(Dose.Info$amt))
          {   ReturnVal <- tkmessageBox(title = "Warning",
                              message = "need DOSE information", icon = "error", type = "ok")
          } else
          {		model<<-"1-comp First Order absorption-Linear"
              param.list<<-c("logV","logKa","logKe")
              One.FirstOrder.linear(Comp1.FO.L,I.Comp1.FO.L.optim,a,b,3,model)     
 			    }
   		    } else if(comp=="one"&admin=="IV bolus"&elim=="linear")
      ## One-comp / IV bolus
   		    {  	if(!is.null(Dose.Info$amt))
              {  ReturnVal <- tkmessageBox(title = "Warning",
                               message = "need DOSE information", icon = "error", type = "ok")
              } else
              {  model<<-"1-comp IV bolus -Linear"
                 param.list<<-c("logV","logKe")
      			     One.IVbolus.linear(Comp1.IVbolus.L,I.Comp1.IVbolus.L.optim,a,b,2,model)    
			        }
   		    } else if(comp=="one"&admin=="IV infusion"&elim=="linear")
     ## One-comp / IV infusion
   	    	{  	if(is.null(Dose.Info$Tinf))
              {   ReturnVal <- tkmessageBox(title = "Warning",
                       message = "need AMT and RATE information", icon = "error", type = "ok")
              } else
              {   model<<-"1-comp IV infusion -Linear"
                  param.list<<-c("logV","logKe")
      			      One.IVinf.linear(Comp1.IVinf.L,I.Comp1.IVinf.L.optim,a,b,2,model)    
              }
   		    } else if(comp=="one"&admin=="Zero Order absorption"&elim=="linear")
     ## One-comp / zero order absorption
   	  	  {  	if(!is.null(Dose.Info$amt))
              {   ReturnVal <- tkmessageBox(title = "Warning",
                               message = "need DOSE information", icon = "error", type = "ok")
              } else
              {		model<<-"1-comp Zero Order absorption -Linear"
                  param.list<<-c("logV","logKe")
 			            tt<-tktoplevel()
			            tktitle(tt)<-"Enter Tk0 parameter"
			            done <- tclVar(0)
			            textEntryWidget <- tkentry(tt, width = paste(50),textvariable = done)
			            tkgrid(tklabel(tt, text = "       "))
			            tkgrid(tklabel(tt, text = "Tk0 :"), textEntryWidget)
			            tkgrid(tklabel(tt, text = "       "))
			            OK.but <- tkbutton(tt, text = "  OK  ",
			                               command = function() 
			                                {  Tk0<<-as.numeric(tclvalue(done));
                                         tkdestroy(tt)
                                         One.ZeroOrder.linear(Comp1.ZO.L,
                                                     I.Comp1.ZO.L.optim,a,b,2,model)  
       					                      })
       					  Cancel.but <- tkbutton(tt, text = "Cancel", command = function() tkdestroy(tt))
    					    tkgrid(OK.but, Cancel.but)
			        }
   		    } else if(comp=="two"&admin=="First Order absorption"&elim=="linear")
     ## Two-comp / First Order Abroption
   		    {	 if(!is.null(Dose.Info$amt))
             {   ReturnVal <- tkmessageBox(title = "Warning",
                               message = "need DOSE information", icon = "error", type = "ok")
             } else
             {  	model<<-"2-comp First Order absorption-Linear"
                  param.list<<-c("logV","logKa","logKe","logK12","logK21")
                  Two.FirstOrder.linear(Comp2.FO.L,I.Comp2.FO.L.optim,a,b,5,model)                                              
			       }
   		    } else if(comp=="two"&admin=="IV bolus"&elim=="linear")
    ## Two-comp / IV bolus
   		    {   if(!is.null(Dose.Info$amt))
              {   ReturnVal <- tkmessageBox(title = "Warning",
                               message = "need DOSE information", icon = "error", type = "ok")
              } else
              {   model<<-"2-comp IV bolus -Linear"
                  param.list<<-c("logV","logKe","logK12","logK21")
      			      Two.IVbolus.linear(Comp2.IVbolus.L,I.Comp2.IVbolus.L.optim,a,b,4,model)    
 			        }
   		    } else if(comp=="two"&admin=="IV infusion"&elim=="linear")
    ## Two-comp / IV infusion
   		    {  	if(is.null(Dose.Info$Tinf))
              {    ReturnVal <- tkmessageBox(title = "Warning",
                       message = "need AMT and RATE information", icon = "error", type = "ok")
              } else
              {    model<<-"2-comp IV infusion -Linear"
                   param.list<<-c("logV","logKe","logK12","logK21")
      			       Two.IVinf.linear(Comp2.IVinf.L,I.Comp2.IVinf.L.optim,a,b,4,model)    
              }
   		    } else if(comp=="two"&admin=="Zero Order absorption"&elim=="linear")
    ## Two-comp / zero order absorption
   		    {  	if(!is.null(Dose.Info$amt))
              {    ReturnVal <- tkmessageBox(title = "Warning",
                               message = "need DOSE information", icon = "error", type = "ok")
              } else
              { 	 model<<-"2-comp Zero Order absorption -Linear"
                   param.list<<-c("logV","logKe","logK12","logK21")
   			           tt<-tktoplevel()
			             tktitle(tt)<-"Enter Tk0 parameter"
			             done <- tclVar(0)
			             textEntryWidget <- tkentry(tt, width = paste(50),textvariable = done)
			             tkgrid(tklabel(tt, text = "       "))
			             tkgrid(tklabel(tt, text = "Tk0 :"), textEntryWidget)
			             tkgrid(tklabel(tt, text = "       "))
			             OK.but <- tkbutton(tt, text = "  OK  ",		
			                       command = function() 
			                       {   Tk0<<-as.numeric(tclvalue(done));
                                 tkdestroy(tt)
                                 Two.ZeroOrder.linear(Comp2.ZO.L,I.Comp2.ZO.L.optim,a,b,4,model)  
       					             })
					         Cancel.but <- tkbutton(tt, text = "Cancel",command = function() tkdestroy(tt))
					         tkgrid(OK.but, Cancel.but)
			        }
   		    } else if(comp=="three"&admin=="First Order absorption"&elim=="linear")
    ## Three-comp / First Order Abroption
   		    {	  if(!is.null(Dose.Info$amt))
              {    ReturnVal <- tkmessageBox(title = "Warning",
                               message = "need DOSE information", icon = "error", type = "ok")
              } else
              { 		model<<-"3-comp First Order absorption-Linear"
                    param.list<<-c("logV","logKa","logKe","logK12","logK21","logK13","logK31")
                    Three.FirstOrder.linear(Comp3.FO.L,I.Comp3.FO.L.optim,a,b,7,model)                                                
			        }
   		    } else if(comp=="three"&admin=="IV bolus"&elim=="linear")
    ## Three-comp / IV bolus
   		    { 	if(!is.null(Dose.Info$amt))
              {    ReturnVal <- tkmessageBox(title = "Warning",
                               message = "need DOSE information", icon = "error", type = "ok")
              } else
              {    model<<-"3-comp IV bolus -Linear"
                   param.list<<-c("logV","logKe","logK12","logK21","logK13","logK31")
      			       Three.IVbolus.linear(Comp3.IVbolus.L,I.Comp3.IVbolus.L.optim,a,b,6,model)    
			        }
   		    } else if(comp=="three"&admin=="IV infusion"&elim=="linear")
    ## Three-comp / IV infusion
   		    {  	if(is.null(Dose.Info$Tinf))
              {   ReturnVal <- tkmessageBox(title = "Warning",
                       message = "need AMT and RATE information", icon = "error", type = "ok")
              } else
              {   model<<-"3-comp IV infusion -Linear"
                  param.list<<-c("logV","logKe","logK12","logK21","logK13","logK31")
      			      Three.IVinf.linear(Comp3.IVinf.L,I.Comp3.IVinf.L.optim,a,b,6,model)    
              }
   		    } else if(comp=="three"&admin=="Zero Order absorption"&elim=="linear")
    ## Three-comp / zero order absorption
   		    {  	if(!is.null(Dose.Info$amt))
              {   ReturnVal <- tkmessageBox(title = "Warning",
                               message = "need DOSE information", icon = "error", type = "ok")
              } else
              {  	model<<-"3-comp Zero Order absorption -Linear"
                  param.list<<-c("logV","logKe","logK12","logK21","logK13","logK31")
 			            tt<-tktoplevel()
			            tktitle(tt)<-"Enter Tk0 parameter"
			            done <- tclVar(0)
			            textEntryWidget <- tkentry(tt, width = paste(50),textvariable = done)
			            tkgrid(tklabel(tt, text = "       "))
			            tkgrid(tklabel(tt, text = "Tk0 :"), textEntryWidget)
			            tkgrid(tklabel(tt, text = "       "))
			            OK.but <- tkbutton(tt, text = "  OK  ",		
			                      command = function() 
			                      {   Tk0<<-as.numeric(tclvalue(done));
                                tkdestroy(tt)
                                Three.ZeroOrder.linear(Comp3.ZO.L,I.Comp3.ZO.L.optim,a,b,6,model)  
       					            })
					       Cancel.but <- tkbutton(tt, text = "Cancel",command = function() tkdestroy(tt))
					       tkgrid(OK.but, Cancel.but)
			        }
   		    } else
   		    {    result<<-NULL
   		    }
	   }		

     RUN.id<<-0
  	 OUTPUT.tot<<-list()
	   PKSAEM.win<<-gwindow("PKmodelFinder",width=300,height=300)
   	 menu.list<-list(File=list( 'Open Data File'=list(handler=OpenDataFile,icon="open"),
                            'Open Dose File'=list(handler=OpenDoseFile,icon="open"),
                            'Quit'=list(handler=function(h,...) dispose(NONMEM.win),icon="quit")),    
                    'EDA'=list('DV vs TIME by ID '=list(handler=ID.plot)),
                    'RUN'=list('Run SAEM'=list(handler=runPK)),
                    'Model Check Plot'=list('DV vs PRED'=list(handler=DVvsPRED.plot),
                                            'DV vs RES'=list(handler=DVvsRES.plot),
                                            'TIME vs RES'=list(handler=TIMEvsRES.plot),
                                            'TIME vs DV and PRED'=list(handler=TIMEvsDVandPRED.plot)),
                    'Covariates'=list('EBE vs COV'=list(handler=EBEvsCOV.plot)))
     gmenu(menu.list,cont=PKSAEM.win)
	   level1<<-gradio(c("one","two","three"))
	   level2<<-gradio(c("IV bolus","IV infusion","First Order absorption","Zero Order absorption"))
	   level3<<-gradio(c("linear","Michaelis Menten"))
	   level4<<-gradio(c("combined error","constant error"))
	   Group1<-ggroup(cont=PKSAEM.win)
     group<-ggroup(horizontal=FALSE,cont=Group1)
     tmp1<-gframe("Compartment",container=group)
     add(tmp1,level1)
     tmp2<-gframe("Administration",container=group)
     add(tmp2,level2)

	   Group2<-ggroup(cont=PKSAEM.win)
     group<-ggroup(horizontal=FALSE,cont=Group2)
     tmp3<-gframe("Elimination",container=group)
     add(tmp3,level3)
     tmp4<-gframe("Error",container=group)
     add(tmp4,level4)
     KWidget <<- gedit("10")
     tmp5<-gframe("# of replication (K)",container=group)
     add(tmp5,KWidget)
     mWidget <<- gedit("500")
     tmp6<-gframe("# of sample (m)",container=group)
     add(tmp6,mWidget)
}

