# Definition of functions .....................................................................

# RUN_FSC2
# run fastsimcoal2
# INPUT
#   input_settings : list with all the information required to run fastsimcoal2
run_fsc2 <- function(input_settings) {
  
  # check that all the required info was given as input
  if(is.null(input_settings$pathToFsc)) {
    stop("Need to specify path to fastsimcoal2 executable file.")
  }
  if(is.null(input_settings$pathTo_InputFile)) {
    stop("Need to specify path to input files.")
  }
  if(is.null(input_settings$TPL_EST_file_tag)) {
    stop("Need to specify tag for TPL, EST and OBS SFS files.")
  }
  if(is.null(input_settings$n_coalsims)) {
    stop("Need to specify number of coalescent simulations.")
  }
  if(is.null(input_settings$n_cycles)) {
    stop("Need to specify number of cycles.")
  }
  
  # copy fsc2 to the working directory
  system(paste("cp ", input_settings$pathToFsc, "fsc26 ", input_settings$pathTo_InputFile,";
               chmod +x ", input_settings$pathTo_InputFile, "fsc26", sep=""))
  
  # fastsimcoal2 default options
  fsc2 <- paste("cd ", input_settings$pathTo_InputFile, ";",
                "./fsc26 ", 
                " -t ", input_settings$TPL_EST_file_tag, ".tpl", 
                " -e ", input_settings$TPL_EST_file_tag, ".est",
                " -L ", input_settings$n_cycles,
                " -n ", format(input_settings$n_coalsims, scientific = F), 
                " -d -M -q -C1 -c2 -B2;", sep="")
  
  
  system(fsc2)
}


# GET_FSC2_COMMANDLINE
# get fastsimcoal2 command line
# INPUT
#   input_settings : list with all the information required to run fastsimcoal2
# RETURN
#   string with the command line
get_fsc2_commandline <- function(input_settings) {
  
  # check that all the required info was given as input
  if(is.null(input_settings$pathToFsc)) {
    stop("Need to specify path to fastsimcoal2 executable file.")
  }
  if(is.null(input_settings$pathTo_InputFile)) {
    stop("Need to specify path to input files.")
  }
  if(is.null(input_settings$TPL_EST_file_tag)) {
    stop("Need to specify tag for TPL, EST and OBS SFS files.")
  }
  if(is.null(input_settings$n_coalsims)) {
    stop("Need to specify number of coalescent simulations.")
  }
  if(is.null(input_settings$n_cycles)) {
    stop("Need to specify number of cycles.")
  }
  
  # fastsimcoal2 default options
  fsc2 <- paste("./fsc26 ", 
                " -t ", input_settings$TPL_EST_file_tag, ".tpl", 
                " -e ", input_settings$TPL_EST_file_tag, ".est",
                " -L ", input_settings$n_cycles,
                " -n ", format(input_settings$n_coalsims, scientific = F), 
                " -d -M -q -C1 -c2 -B2;", sep="")
  
  fsc2
}



# PLOT_FITSFS_1D
# plot the 1d sfs fit to the data
#   obs.sfs : observed SFS
#   exp.sfs : expected SFS
#   ylims   : limits for the y axis
#   pop.name : pop.name
plot_fitSFS_1d <- function(obs.sfs, exp.sfs, ylims, pop.name) {
  par(mfrow=c(1,2), oma=c(2,2,5,2))
  mymat <- matrix(c(as.numeric(obs.sfs), as.numeric(exp.sfs)*sum(obs.sfs)), ncol=2)
  barplot(t(mymat), beside = TRUE, names.arg = 0:(length(obs.sfs)-1), legend.text =c("obs","exp"), ylab="count SNPs", main="Fit SFS")
  # plot(0:(length(obs.sfs)-1), obs.sfs, type="s", ylim=ylims, lwd=2, xlab=pop.name, ylab="SNP counts")
  # lines(0:(length(obs.sfs)-1), exp.sfs*sum(obs.sfs), type="s", col=3, lty=1, lwd=2)
  # legend("topright", c("obs","exp"), col=c(1,3), lty=c(1,1), lwd=2)    
  # plot(0:(length(obs.sfs)-1), log10(obs.sfs), type="s", ylim=log10(ylims), lwd=2, xlab=pop.name, ylab="log10(SNPcounts)")
  barplot(t(log10(mymat)), beside = TRUE, names.arg = 0:(length(obs.sfs)-1), ylab="log10(count SNPs)", main="Fit SFS log10 scale", ylim=log10(c(min(mymat[mymat>0])*0.95,max(mymat)*1.05)), xpd = FALSE)
  # lines(0:(length(obs.sfs)-1), log10(exp.sfs*sum(obs.sfs)), type="s", col=3, lty=1, lwd=2)
  # legend("topright", c("obs","exp"), col=c(1,3), lty=c(1,1), lwd=2)
}


########################
##
## AIC related functions - Laurent
##  (c) Laurent Excoffier Computation of AIC and relative lhoods
##
#######################


log10toln<-function(l10) {
  logn=l10/log10(exp(1))
}
#...............................................................................

AIC<-function(log10L, k) {
  lognL=log10toln(log10L)
  2*k-2*lognL
}
#...............................................................................
relLhood<-function(AICs) {
  minAIC=min(AICs)
  ws=exp(0.5*(minAIC-AICs))
  ws
}

# COMPUTELHOOD
# computes the loglikelihood given an obs.sfs and exp.sfs
# INPUT
#   obs.sfs : numeric vector with the observed SFS
#   exp.sfs : numeric vector with the expected SFS
# RETURN
#   log likelihood computed as SUM m_i log10(p_i), 
#   where m_i is the ith entry of obs.sfs and p_i is the ith entry of exp.sfs
# NOTE: for entries where m_i > 0 and p_i=0, we replace p_i by a small value (penalty)
computelhood <- function(obs.sfs, exp.sfs) {
  
  lhood <- 0
  
  # remove the first and last entries
  obs.sfs <- obs.sfs[-c(1,length(obs.sfs))]
  exp.sfs <- exp.sfs[-c(1,length(exp.sfs))]
  
  # Get the valid entries, i.e. entries where obs.SFS > 0
  eval <- which(obs.sfs > 0)
  
  # Calculate expected SFS with the penaltie for entries where obs.SFS > 0 and exp.SFS == 0
  if(sum(exp.sfs[eval]==0) > 0) {
    # Settings (penalty for exo SFS entries with zero)
    penalty <- 1e-10
    minpenalty <- 1e-8
    
    penalty <- min(exp.sfs[exp.sfs>0])/100
    if(penalty > minpenalty) {
      penalty <- minpenalty
    } 
    
    # Get the entries which are zero at the obs SFS to have the penalty
    tmp.exp <- exp.sfs[eval]  # note that the length of tmp.exp is length(eval)
    tmp.exp[tmp.exp==0] <- penalty 
    exp.sfs[eval] <- tmp.exp
  }
  
  # check that the sum of exp.sfs is 1
  exp.sfs <- exp.sfs/sum(exp.sfs)
  
  # compute the likelihood
  if(sum(exp.sfs[eval]==0) > 0) {
    print("ERROR: still entries with exp.sfs==0!!!!")
  } else {
    lhood <- sum(obs.sfs[eval]*log10(exp.sfs[eval]))    
  }
  
  lhood
}


# Function to add error bars - copied from http://monkeysuncle.stanford.edu/?p=485
error.bar <- function(x, upper, lower=upper, length=0.1,...){
  if(length(x) !=length(lower) | length(lower) != length(upper)) {
    stop("vectors must be same length")
  }    
  arrows(x,upper, x, lower, angle=90, code=3, length=length, ...)
}


#--- Function to remove separators within a string
removeTrailingSep=function(string, sep) {
  temp=strsplit(string, split=sep)
  temp2=temp[[1]][nchar(temp[[1]])>0]
  cleanStr=temp2[1]
  if (length(temp2)>1) {
    for (i in 2:length(temp2)) {
      cleanStr=paste(cleanStr, temp2[i], sep=sep)
    }
  }
  cleanStr
}

#--- Replace Keep by -9999
replaceKeep=function(string) {
  if (grepl("keep", string)) {
    return(gsub("keep", "-9999", string))
  }
  return(string)
}

#--- Reading numbers on separate lines -----
getNumbers=function(start, parFile, numSamples) {  
  for (i in 1:numSamples) {
    curnum=as.numeric(unlist(strsplit(parFile[start+i], split=separator))[1])
    if (i==1) {
      num=curnum 
    } else {
      num=c(num, curnum)
    }
  }
  num
}

readSampleSizesTimesAndInbreedingLevel=function(start, parFile, numSamples) {
  for (i in 1:numSamples) {
    curLine=unlist(strsplit(parFile[start+i], split=separator))
    curSampSize=as.numeric(curLine[1])
    curSampTime=0
    curInbreeding=0
    if (length(curLine)>1) curSampTime=as.numeric(curLine[2])
    if (length(curLine)>2) curInbreeding=as.numeric(curLine[3])
    if (i==1) {
      sampSize=curSampSize
      sampTime=curSampTime
      inbreeding=curInbreeding
    } else {
      sampSize=c(sampSize,curSampSize)
      sampTime=c(sampTime,curSampTime)
      inbreeding=c(inbreeding,curInbreeding)
    }
  }
  list(ss=sampSize, st=sampTime, inb=inbreeding)
}

#--- Read migration matrix
readMigMat=function(start, parFile, numSamples) {
  for (i in 1:numSamples) {
    tmpseparator=separator
    if(is.na(as.numeric(unlist(strsplit(parFile[start+i], split=separator)))[1])) {
      tmpseparator=separator2
    }
    curmigs=as.numeric(unlist(strsplit(parFile[start+i], split=tmpseparator)))
    if (i==1) {
      migs=curmigs 
    } else {
      migs=rbind(migs, curmigs)
    }
  }
  rownames(migs)=1:numSamples
  migs 
}

# SIMULATE_EXPECTEDSFS
# simulate the expected SFS according to model parameters
# INPUT
#   input_settings : list with all the information required to run fastsimcoal2
simulate_expectedSFS_1d <- function(input_settings, obs_file_end) {
  obsfilename <- paste(input_settings$TPL_EST_file_tag, obs_file_end, sep="")
  system(paste("cd ", input_settings$pathTo_InputFile, ";", "cp ", obsfilename, " ", input_settings$TPL_EST_file_tag, ";", sep=""))
  
  system(paste("cd ", input_settings$pathTo_InputFile, input_settings$TPL_EST_file_tag, ";
                mv ", obsfilename, " ", input_settings$TPL_EST_file_tag, "_maxL", obs_file_end, ";",  sep=""))
  
  system(paste("cd ", input_settings$pathTo_InputFile, input_settings$TPL_EST_file_tag,  ";", 
               "../fsc25221 -i ", input_settings$TPL_EST_file_tag, "_maxL.par -n ", format(input_settings$n_coalsims, scientific=FALSE)," -d -q -c4 -B4;
              rm ", input_settings$TPL_EST_file_tag, "_maxL", obs_file_end, ";", sep=""))  
}


# DERIVED2MAF
# transforms 2D derived SFS into 2D MAF SFS
# INPUT:
#   derived_sfs2d : 2d matrix with the 2d derived SFS
# RETURN
#   maf_sfs2d : 2d matrix with the 2d MAF SFS
derived2maf <- function(derived_sfs2d) {
  
  n1 <- nrow(derived_sfs2d)
  n2 <- ncol(derived_sfs2d)
  
  maf_2dsfs <- matrix(0, nrow=n1, ncol=n2)
  
  threshold_freq <- 0.5*(n1+n2-2)
  
  for(i in 0:(n1-1)) {
    for(j in 0:(n2-1)) {
      if(i+j < threshold_freq) {
        maf_2dsfs[i+1, j+1] <-  maf_2dsfs[i+1, j+1]+derived_sfs2d[i+1, j+1]
      } 
      else if(i+j == threshold_freq) {
        maf_2dsfs[i+1, j+1] <- maf_2dsfs[i+1, j+1]+(0.5*derived_sfs2d[i+1, j+1]+0.5*derived_sfs2d[n1-i,n2-j])           
      }
      else {
        maf_2dsfs[n1-i, n2-j] <- maf_2dsfs[n1-i, n2-j]+derived_sfs2d[i+1, j+1]
      }
    }     
  }
  maf_2dsfs
}


# Plot 2D SFS for observed and expected SFS
plot2dSFS <- function(obsSFS,expSFS,xtag,ytag, minentry) {
  layout(matrix(1:3, nrow = 1), widths = c(0.4,0.4,0.2))
  
  # transform data.frames into matrices
  obsSFS <- as.matrix(obsSFS)
  expSFS <- as.matrix(expSFS)
  
  library(RColorBrewer)
  nclasses <- c(0,colorRampPalette(brewer.pal(8,"OrRd"))(15))
  
  obsSFS[1,1] <- 0
  expSFS <- expSFS/sum(expSFS)
  expSFS <- round(expSFS*sum(obsSFS))
  
  breaksplot <- c(0, log10(minentry),seq(log10(minentry),max(log10(obsSFS)),length.out = length(nclasses)-1))
  image(1:nrow(obsSFS), 1:ncol(obsSFS), log10(obsSFS), col=nclasses, breaks=breaksplot, xlab=xtag, ylab=ytag, main="Observed 2D-SFS", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  image(1:nrow(expSFS), 1:ncol(expSFS), log10(expSFS), col=nclasses, breaks=breaksplot, xlab=xtag, ylab=ytag, main="Expected 2D-SFS", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  image.scale(log10(obsSFS), zlim=range(breaksplot), col = nclasses, breaks=breaksplot, horiz=FALSE, ylab="log10(SFS counts)", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  
}

# Plot fit 2D SFS for the relative SFS
plot_relDiff2dSFS <- function(obsSFS,expSFS,xtag,ytag, minentry) {
  layout(matrix(1:2, nrow = 1), widths = c(0.7,0.3))
  
  library(RColorBrewer)
  
  obsSFS[1,1] <- 0
  expSFS <- expSFS*sum(obsSFS)
  
  eval <- obsSFS > minentry
  diffSFS <- matrix(0,nrow=nrow(obsSFS), ncol=ncol(obsSFS))
  diffSFS[eval] <- (obsSFS[eval]-expSFS[eval])/obsSFS[eval]
  eval <- obsSFS <= minentry
  meandiff <- (sum(obsSFS[eval], na.rm=T)-sum(expSFS[eval], na.rm=T))/sum(obsSFS[eval], na.rm=T)
  diffSFS[eval] <- meandiff/sum(eval)
  
  nclasses <- rev(c(colorRampPalette(brewer.pal(8,"BrBG"))(15)))
  nclasses[8] <- 0
  breaksplot <- c(seq(-1,-0.1,length.out = (length(nclasses))/2), seq(0.1,1,length.out = (length(nclasses))/2))
  image(1:nrow(expSFS), 1:ncol(expSFS), diffSFS, col=nclasses, breaks=breaksplot, xlab=xtag, ylab=ytag, main="Relative difference obs.-exp.", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  image.scale(diffSFS, zlim=range(breaksplot), col = nclasses, breaks=breaksplot, horiz=FALSE, ylab="(obs. SFS - exp. SFS)/obs. SFS", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  
}

#This function creates a color scale for use with the image()
#function. Input parameters should be consistent with those
#used in the corresponding image plot. The "horiz" argument
#defines whether the scale is horizonal(=TRUE) or vertical(=FALSE).
image.scale <- function(z, zlim, col = rainbow(12), breaks, horiz=TRUE, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){ylim<-c(0,1); xlim<-range(breaks)}
  if(!horiz){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}




