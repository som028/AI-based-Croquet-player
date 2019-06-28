#Genetic Algorithm

#population defination
pop<-function(pop_size=50,Lx=0,Ux=1500,Ly=0,Uy=3000,CHro_size=100)
{
  newpop=chromosome(Lx,Ux,Ly,Uy,CHro_size) #creating chromosome of (x,y) coordinate
  for(i in c(1:(pop_size-1))) {
      newpop=cbind(newpop,chromosome(Lx,Ux,Ly,Uy,CHro_size)) }
  return(newpop)
}

chromosome<-function(Lx=0,Ux=1500,Ly=0,Uy=3000,Chro_size=100) #randomly creating chromosome within boundary limit
{
  Ch=0
  #Sampling of chromosome length from the solution space
  Ch=data.frame("x"=sample(seq(Lx+400,Ux-400,4),Chro_size),"y"=sample(seq(Ly+200,Uy-600,4),Chro_size))
  return(Ch)
}


feild<-function(Lx=0,Ux=1500,Ly=0,Uy=3000)
{
  x=(Ux-Lx)/2
  y=(Uy-Ly)*2/30
  po=data.frame("x"=x,"y"=y) #start point
  
  #forward feild 
  po=rbind(po,c(x,y*2))
  po=rbind(po,c(x,y*3))
  po=rbind(po,c(Ux-400,y*5))
  po=rbind(po,c(x,y*7))
  po=rbind(po,c(Ux-400,y*9))
  po=rbind(po,c(x,y*11))
  po=rbind(po,c(x,y*12))
  
  #backward feild coordinate
  po=rbind(po,c(x,y*12))
  po=rbind(po,c(x,y*11))
  po=rbind(po,c(Lx+400,y*9))
  po=rbind(po,c(x,y*7))
  po=rbind(po,c(Lx+400,y*5))
  po=rbind(po,c(x,y*3))
  po=rbind(po,c(x,y*2))
  # plot(po,xlim=c(Lx,Ux))
  # for(i in c(2:nrow(po)))
  #   segments(po[i,1]-30,po[i,2],po[i,1]+30,po[i,2],lwd=4,col="blue")
  return(po)
}

line<-function(pol,Ini,Nex,l)
{
  deltay=(Nex[,2]-Ini[,2]) #detla x
  deltax=Nex[,1]-Ini[,1] #delta y
  if(abs(deltax)<0.01) #checking if delta x is greater then zero non-infinite slope
    deltax=0.01
  m=deltay/deltax #Slope m
  c=(Ini[,2]-m*Ini[,1]) #Intercept c
  d1=(pol[,2]-l-m*(pol[,1]-l)-c)/(sqrt(1+m**2)) #Distance from first pole of wicket
  d2=(pol[,2]+l-m*(pol[,1]+l)-c)/(sqrt(1+m**2)) #
  d1=abs(d1)
  d2=abs(d2)
  #print(d1)
  if(d1<31 | d2<31)
    return(1)
  else
    return(0)
}

reward<-function(popu,Lx=0,Ux=1500,Ly=0,Uy=3000)
{
  pass=rew=0
  rew[1:(ncol(popu)/2)]=0
  pass[1:(ncol(popu)/2)]=0
  pol_arr=feild()
  
  
  for(index in c(1:(ncol(popu)/2)))
  {
    initial = pol_arr[1,]
    curr=2
    j=index*2-1
    for(i in c(1:nrow(popu)))
    {
      bx=popu[i,j]
      by=popu[i,j+1]
      if(by<Ly | by>Uy | bx<Lx | bx>Ux) #outside boundary penalty
        rew[index]=rew[index]-0.5
      
      px=pol_arr[curr,1]
      py=pol_arr[curr,2] #y-coordinate of pole
      
      # if(index>16)
      # {
      #   print(pol_arr[curr,])
      #   print(initial)
      #   print(popu[i,j:(j+1)])
      # }
      
      dis=line(pol=pol_arr[curr,],Ini=initial,Nex=popu[i,j:(j+1)],30) #check the range for ball hit
      
      initial=popu[i,j:(j+1)] #initial point of ball
      
      if(curr<8)
      {
        if(dis==1)
        {
          if(by>py)
          {
            rew[index]=rew[index]+1
            curr=curr+1
          }
        }
        else
        {
          if(by>py)
          {
            rew[index]=rew[index]-0.5
            #curr=curr+1
          }
        }
      } else
      {
        if(dis==1)
        {
          if(by<py)
          {
            rew[index]=rew[index]+1
            curr=curr+1
          }
        }
        else
        {
          if(by<py)
          {
            rew[index]=rew[index]-0.5
            #curr=curr+1
          }
        }
      }
      if(curr>=nrow(pol_arr)) #reach the end
      {
        rew[index]=rew[index]*100/i #proportional reward on finishing
        pass[index]=i
        break
      }
    }
    #print(curr)
    if(curr<nrow(pol_arr)) #did not reach the end
    {
      rew[index]=rew[index]-0.5
      pass[index]=i
    } 
  }
  #rew=(rew/(sum(abs(rew))))*100 #scaling of reward
  return(list("re"=rew,"hit"=pass))
}

crosover<-function(p1,p2,Crosrate=0.7) #preparing new chromosome by exchanging the parents chromosomes
{
  Cr=sample(c(1:nrow(p1)),Crosrate*100)
  a=p1
  p1[Cr,]=p2[Cr,]
  p2[Cr,]=a[Cr,]
  child=cbind(p1,p2) 
  return(child)
}

mutate<-function(p,Mutrate=0.05,Lx=0,Ux=1500,Ly=0,Uy=3000) #mutation by exchanging the order of chromosomes
{
  lx=Ux-Lx #mutation limit for x
  ly=Uy-Ly #mutation limit for x
  mv=round(Mutrate*100) #nuumber of values to mutate
  muix=sample(c(1:nrow(p)),mv) #randomly select the index for mutation for x
  muiy=sample(c(1:nrow(p)),mv) #randomly select the index for mutation for x
  mux=sample(c(-lx:lx)*Mutrate,mv) #Choosing random values to mutate x between limit
  muy=sample(c(-ly:ly)*Mutrate,mv) #Choosing random values to mutate y between limit
  p$x[muix]<-p$x[muix]+mux #mutation of x coordinate
  p$y[muiy]<-p$y[muiy]+muy #mutation of y coordinate
  
  #if coordinateKeeping coordinates are outside the boundary putting a limit on them
  p$x[p$x>Ux]=Ux+10
  p$x[p$x<Lx]=Lx-10
  p$y[p$y>Uy]=Uy+10
  p$y[p$y>Uy]=Ly-10
  
  #returning the muatated chromosome
  return(p)
}

selectionpool<-function(tourpop,k=3)
{
  #for defining fitness of each individual
  newpop=data.frame(c(1:nrow(tourpop)))
  i=index=k2=p=0
  p[1:k]<-0
  k3=1
  Popsize=ncol(tourpop)/2
  for (j in seq(1,Popsize,1))
  {
    k2=0
    p[1:k]<-0
    p<-p[1:k]
    index=sample(c(1:Popsize),k)
    i=index*2-1
    for(k1 in i)
    {
      k2=k2+1
      p[k2]=reward(tourpop[,k1:(k+1)])$re
      if(k2>1)
      {
        if(p[k2]>p[k3]) {
          k3=k2 }
      }
    }
    #print(p)
    newpop=cbind(newpop,tourpop[,i[k3]:(i[k3]+1)])
  }
  newpop=newpop[,2:ncol(newpop)]
  return(newpop)
}

#New generation
makegen<-function(tourpop, crosrate=0.5, mutrate=0.05) #Crossover and mutation occurs here
{
  newpop=data.frame(c(1:nrow(tourpop))) #Intialize the population
  gen=i=0
  for(index in seq(1,(ncol(tourpop)/2),2)) #For loop 
  {
    i=index*2-1 #position in popultaion
    P1=tourpop[,i:(i+1)] #Parent 1
    P2=tourpop[,(i+2):(i+3)] #Parent 2
    gen=crosover(p1=P1,p2=P2) #crossover of two parents
    gen[,1:2]=mutate(gen[,1:2]) #mutaion of child 1
    gen[,3:4]=mutate(gen[,3:4]) #mutaion of child 2
    newpop=cbind(newpop,gen)
  }
  newpop=round(newpop[,2:ncol(newpop)],2)
  
  #Elitism  
  e=round(ncol(tourpop)/20) #selecting 10% elite individuals
  fit=data.frame("index"=c(1:(ncol(tourpop)/2)),"fi"=reward(tourpop)$re) #Index and reward
  elite=fit[order(-fit$fi),] 
  elite=elite[1:e,]
  elite$index=elite$index*2-1
  elind=c(elite$index,elite$index+1)
  elind=elind[order(elind)] 
  elipop=tourpop[,elind]
  newpop[,1:ncol(elipop)]=elipop #Replace the population with elite population
  
  return(newpop)
}

GeneticEvolution<-function(lx=0,ux=1500,ly=0,uy=3000, Chrosize=140,Popsize=120, Gen=30, Crosrate=0.7, Mutrate=0.05,K=3)
{
  firsc=-100
  
  popu = pop(pop_size=Popsize,Lx=lx,Ux=ux,Ly=ly,Uy=uy,CHro_size=Chrosize)  #population initialization
  cl=colnames(popu)
  fit=reward(popu)$re                    #fitnes of each individual
  hits=reward(popu)$hit
  maxfit=data.frame("Points"=max(fit),"Hits"=min(hits)) #maximum fitnes of initial population
  for(j in 1:Gen)
  {
    popu=selectionpool(popu,k=K)    #New selection pool
    print(j)
    popu=makegen(popu, crosrate=Crosrate, mutrate=Mutrate) #New Generation
    
    ree=reward(popu) #calculating points and hits of population
    rr=ree$re #calacluating fitness of each individual
    hh=ree$hit #calacluating hitness of each individual
    
    
    mf=max(rr) #calacluating maximum fitness of population
    mh=min(hh) #calacluating minimum hitness of population
    
    if(mf>firsc) #comparing the fitness 
    {
      firsc=mf #Replacing the variable with higher fitness value
      fitgen=j+1 #Storing the generation number of maximum fitness
      fitpos=which.max(rr) #reteriving the position of maximum fitness
      #print(fitpos)
      fh=(hh[fitpos]) #Hits taken by player to reach that fitness value
      fitpos=fitpos*2-1
      
      bestfitplayer<-popu[1:fh,fitpos:(fitpos+1)] #Storing the x and y coordinate of player with best fitness
      
    }
    
    maxfit=rbind(maxfit,c(mf,mh)) #maximum fitness for each generation
    row.names(popu)<-c(1:Chrosize) 
    colnames(popu)<-cl
    
  }
  

  par(mfrow=c(2,2))
  
  write.csv(bestfitplayer,file="bestfitplayer.csv")
  print(bestfitplayer)
  
  #Plotting max fitness at each generaiton
  fitness=data.frame("Generation"=seq(1,Gen+1,1),"Reward"=maxfit$Points)
  plot(fitness) 
  lines(with(fitness,smooth.spline(Generation, Reward, spar=0.35)))
  abline(v=fitgen,col='red')
  
  #Plotting min hits at each generation
  hitness=data.frame("Generation"=seq(1,Gen+1,1),"Hits"=maxfit$Hits)
  plot(hitness)
  lines(with(hitness,smooth.spline(Generation, Hits, spar=0.35)))
  abline(v=fitgen,col='red')
  
  #Plotting the feild and path
  po=feild()
  plot(po,xlim=c(Lx,Ux))
  for(i in c(2:nrow(po)))
     segments(po[i,1]-30,po[i,2],po[i,1]+30,po[i,2],lwd=4,col="blue")
  par(new=TRUE)
  plot(bestfitplayer)
  lines(bestfitplayer)
  
  #return(bestfitplayer)
}

GeneticEvolution() #Running the function


