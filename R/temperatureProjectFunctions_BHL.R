# sample jags with chunks

sampleJagsInChunks <- function(cl,fileOutName,bugsName,data,inits,nAdapt,params,nIterChunk,nThin){
  
  ( beforeJags <- Sys.time() )
  print( beforeJags )
  print( 'Adapting, no progress bar available when running chains in parallel' )
  beforeAdapt <- Sys.time()
  ###################################
  # adaptation
  ####################################
  
  adapt <- clusterEvalQ(cl, {
    require(rjags)
    
    jm <- jags.model(
      file = bugsName,#paste0(baseDir,bugsName),         
      data = data,
      inits = inits,
      n.chains = 1,
      n.adapt = nAdapt,             
    )
    return( jm )
  })
  
  afterAdapt <- Sys.time()
  print(paste('adaptation time:',afterAdapt - beforeAdapt,sep=" "))
  adaptTime = afterAdapt - beforeAdapt
  
  ####################################
  # sampling
  ####################################
  
  # run 1 iter just to get structure of outP
  outP <- clusterEvalQ(cl, {
    
    require(rjags)
    load.module("dic")
    
    js <- jags.samples(  ##############coda.samples( #
      model = jm,
      variable.names = params,
      n.iter = 1,
      thin = nThin,
      progress.bar = 'text'
    ) 
    return( js )  ###############as.mcmc
  })
  
  
  print( paste( Sys.time(), Sys.time() - beforeJags, sep=" " ) )
  
  ####################################
  # run separate chunks of iterations, nIterChunk at a time
  # each chunk gets appended to outP and saved as 'fileOutName'
  ####################################
  sampleTime <- array( NA,nSave )  
  
  for( i in 1:nSave ){
    
    holdTime <- Sys.time() 
    print( paste( 'Before js, i = ',i,'; Start time = ',holdTime, sep="" ) ) 
    
    outP[[i]] <- clusterEvalQ(cl, {
      
      require(rjags)
      load.module("dic")    
      
      js <- jags.samples( ##############coda.samples( #
        model = jm,
        variable.names = params,
        n.iter = nIterChunk,
        thin = nThin,
        progress.bar = 'text'
      ) 
      
      return( (js) )  ##############as.mcmc
    })
    
    print( paste( 'After js, Loop = ',i,' out of ',nSave, ', Saved iters done = ',nIterChunk*i, sep="") )
    sampleTime[i] = Sys.time() - holdTime
    print( paste( 'Chunk end time =',Sys.time(), '; Chunk run time =', round( Sys.time() - holdTime, 2 ), sep=" " ) )
    print( getwd() )  
    print( paste( 'mean sampleTime = ', round(mean(sampleTime, na.rm=T),2), '; adaptTime = ', round( adaptTime,2 ), sep='' ) )
    
    print( "##########" )
    
    save(d, outP, i, adaptTime, sampleTime, file = fileOutName)  
  }
  
  
  ( done <- Sys.time() )
  
  print(afterAdapt - beforeAdapt) 
  print(done - beforeJags)
  
  stopCluster(cl)
}
#######################################################


#######################################
# jags samples
#######################################
# append chains together for each chunk
# this doesn't work when just one chunk has been run becasue the trailing lists with just one iter get counted
appendJagsChunks <- function( outP,nThin,nIter,nIterChunk ){
  
  nChunks <- length(outP)
  nChains <- length(outP[[1]]) #3
  nVar <-  length(outP[[1]][[1]]) #
  #nThin <- 5 
  
  #nIter <- 25000                      # total num of iters
  #nIterChunk <- 100 #round(nIter/nSave)   # num of iters per chunk
  nSave <- round(nIter/nIterChunk)                       # num of saves
  nOut <- nIterChunk/nThin 
  
  o2 <- outP #initialize list
  
  for( i in 1:(nChunks-0) ){
    print(c(i))
    for(v in 1:nVar){
      chainIndex <- length(dim(outP[[1]][[1]][[v]])) 
      if(nChains==3) o2[[i]][[v]] <- abind(outP[[i]][[1]][[v]],outP[[i]][[2]][[v]],outP[[i]][[3]][[v]],along=chainIndex)
      if(nChains==4) o2[[i]][[v]] <- abind(outP[[i]][[1]][[v]],outP[[i]][[2]][[v]],outP[[i]][[3]][[v]],outP[[i]][[4]][[v]],along=chainIndex)    
      if(nChains==5) o2[[i]][[v]] <- abind(outP[[i]][[1]][[v]],outP[[i]][[2]][[v]],outP[[i]][[3]][[v]],outP[[i]][[4]][[v]],outP[[i]][[5]][[v]],along=chainIndex)    
    }
  }
  
  o3 <- list()#o2
  varNames <- names(outP[[1]][[1]])
  # append chunks
  for(v in 1:nVar){
    print(v)
    hold <- o2[[1]][[v]]
    for( i in 2:(nChunks-0) ){    
      iterIndex <- length(dim(outP[[1]][[1]][[v]])) - 1
      hold <- abind( hold,o2[[i]][[v]], along=iterIndex )
    }
    #names(hold) <- varNames[v]
    o3[[v]] <- hold
  }    
  out <- o3
  names(out) <- names(outP[[1]][[1]])
  return(out)
}
