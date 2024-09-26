# extra plotting functions

# Shelf to open ocean simulations
shelf_slope_open_sim <- function(data) {
  out <- data.frame()
  sims <- list()
  
  for(i in 1:nrow(data)) {
    
    id <- data$id
    
    my_setup <- setupVertical2(szprod = data[i, 4], # small mesozoopl production
                               lzprod = data[i, 5], # large mesozoopl production
                               bprodin = data[i, 3], # Benthic production
                               depth = as.numeric(data[i, 2]), # Bottom depth
                               photic = 150) # Photic zone depth default
    
    
    # Run FEISTY:
    sim <- simulateFEISTY(p = my_setup) 
    
    # Take average biomass of the last 80 years:
    totBiomass <- sim$totBiomass %>%
      as.data.frame() %>%
      tail(80) %>%
      summarise(across(everything(), mean)) %>%
      rename(smallPel  = totBiomass.smallPel,
             mesoPel   = totBiomass.mesoPel,
             largePel  = totBiomass.largePel,
             midwPred  = totBiomass.midwPred,
             demersals = totBiomass.demersals) %>%
      t() %>%
      as.data.frame() %>%
      rename(totBiomass = V1) %>%
      rownames_to_column(var = "funGroup") %>%
      mutate(id = id[i])
    
    out <- rbind(out, totBiomass)
    sims[[i]] <- sim
    
  }
  
  out <- out %>%
    mutate(totBiomass = case_when(totBiomass == 0 ~ NA,
                                  T ~ totBiomass))
  
  return(list(sims = sims,
              Biomass = out))
}

# network plot function that has a fixed scale for biomass bubbles
plotNetwork2 <- function(sim, scale_bio_manual=F, scale_bio = 5) {
  p <- sim$p
  u <- sim$u
  
  # Number of groups and biomass:
  ngroup <- p$ix[[length(p$ix)]][length(p$ix[[length(p$ix)]])] #resources (4) + fish 
  biomass <- u
  
  #Average of the biomass (mean value of the last 40% time) : 
  Bi <- colMeans(biomass[round(0.6*nrow(biomass), digits = 0):nrow(biomass),]) # 
  
  if (p$setup == "setupBasic" | p$setup == "setupBasic2"){
    
    # Set artificial depths to offset bubbles:
    Av_depth <- c(-1,-1,-4,-4,rep(0, length(p$ix[[1]])),rep(-2, length(p$ix[[2]])),rep(-3, length(p$ix[[3]])))
    
    p$SpId <- p$groupnames[(p$nResources+1):length(p$groupnames)]
    SpId <- c(p$groupnames[1:p$nResources],
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])))
    
    # Specify depth axis text:
    yaxis <- c("", "") # Leave blank because there is no depth in setup vertical
    
    # Set artificial depth:
    p$bottom <- -(min(Av_depth)) + 1
  }  
  
  if (p$setup == "setupVertical" | p$setup == "setupVertical2"){
    
    #Calculate average depth day/night
    
    Av_depth_day <- 1 : p$nStages
    Av_depth_night <- 1 : p$nStages
    for (i in 1:p$nStages) {
      Av_depth_day[i] <- which.max(p$depthDay[ ,i])
      Av_depth_night[i] <- which.max(p$depthNight[ ,i]) 
    }
    
    Av_depth <- -(Av_depth_day + Av_depth_night) / 2
    
    # Change a bit for visualization:
    Av_depth[p$ix[[1]][1]:p$ix[[1]][length(p$ix[[1]])]] <- Av_depth[p$ix[[1]][1]:p$ix[[1]][length(p$ix[[1]])]] + 0.1 * p$bottom
    Av_depth[p$ix[[3]][1]:p$ix[[3]][length(p$ix[[3]])]] <- Av_depth[p$ix[[3]][1]:p$ix[[3]][length(p$ix[[3]])]] - 0.1 * p$bottom
    Av_depth[p$ix[[4]][1]:p$ix[[4]][length(p$ix[[4]])]] <- Av_depth[p$ix[[4]][1]:p$ix[[4]][length(p$ix[[4]])]] - 0.1 * p$bottom
    Av_depth[p$ix[[2]][1]:p$ix[[2]][length(p$ix[[2]])]] <- Av_depth[p$ix[[2]][1]:p$ix[[2]][length(p$ix[[2]])]] - 0.2 * p$bottom
    
    # Set color palette: 
    p$SpId <- p$groupnames[(p$nResources+1):length(p$groupnames)]
    SpId <- c(p$groupnames[1:p$nResources],
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])),
              rep(p$SpId[4], length(p$ix[[4]])),
              rep(p$SpId[5], length(p$ix[[5]])))
    
    # Specify depth axis text:
    yaxis <- c("Epipelagic   ", "     Bottom")
  }
  
  # Marker size depends on biomass following a cubic square transformation
  Max_bio <- ifelse(scale_bio_manual, scale_bio, max(Bi, na.rm = T))
  Msize   <- Bi / Max_bio
  Msize[Msize == 0] <- NA
  if(scale_bio_manual==T) Msize <- 1.2*Msize^(1/2)
  if(scale_bio_manual==F) Msize <- Msize^(1/3)
  Max_size = max(Msize, na.rm = T)
  
  # Create line width: 
  Flux <- FEISTY:::getFeeding(sim)
  Flux <- c(Flux) 
  threshold <- 0.05 
  indx <- which(Flux >= threshold) # takes the x highest values of theta
  
  # Set values of each coordinate (i.e. size and water column position) and put together:
  coord_1 <- data.frame(index = 1:p$nStages^2,
                        mc = rep(p$mc[1:p$nStages], p$nStages), 
                        depth = rep(Av_depth[1:p$nStages], p$nStages), 
                        SpId = rep(SpId, p$nStages),
                        Msize = rep(Msize, p$nStages), 
                        LineWdth = (Flux/max(Flux))^(1/3),
                        Alpha = (Flux/max(Flux))^(1/3))
  
  coord_2 <- data.frame(index = 1:p$nStages^2, # Notice repetition of ys grouped by "each" to change order
                        mc = rep(p$mc[1:p$nStages], each = p$nStages), 
                        depth = rep(Av_depth[1:p$nStages], each = p$nStages), 
                        SpId = rep(SpId, each = p$nStages),
                        Msize = rep(Msize, each = p$nStages),
                        LineWdth = (Flux/max(Flux))^(1/3),
                        Alpha = (Flux/max(Flux))^(1/3))
  
  # Combine in a data frame:
  df <- rbind(coord_1, coord_2)
  df <- df[order(-df$Msize),]   
  df <- subset(df,!(is.na(df$Msize))) 
  df2 <- subset(df,df$index %in% indx)
  df2 <- df2[order(-df2$Msize),]
  
  df$SpId=factor(df$SpId,p$groupnames)
  df2$SpId=factor(df2$SpId,p$groupnames)
  
  # Generate plot:
  plot <- FEISTY:::defaultplot() +
    geom_line(data = df2, aes(x = mc, y = depth, group = index, color = SpId, alpha = Alpha),
              show.legend = F, linewidth = df2$LineWdth) +
    geom_point(data = df, aes(x = mc, y = depth, color = SpId, size = Msize), stroke = 0, shape = 16) +
    scale_color_manual(values = p$my_palette[attr(p$my_palette, "names") %in% df$SpId], 
                       labels = p$my_names[attr(p$my_palette, "names") %in% df$SpId]) +
    scale_radius(limits = c(0, NA), range = c(0, (8* Max_size)))+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_continuous(breaks = seq(0, round(-p$bottom - 1), by = -p$bottom), labels = yaxis) +
    annotation_logticks(sides = "b",linewidth = 0.4,colour = "darkgrey") +
    labs(x ="Weight (g)", y = "", color = "Groups") +
    guides(size = "none",
           color = guide_legend(override.aes = list(size = 5),
                                nrow = 2,
                                byrow = TRUE)) +
    theme(legend.position = "bottom",legend.key = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=8),
          axis.text.y = element_text(angle = 90, hjust = .5, margin = margin(r = 0), size=8),
          axis.ticks.y = element_blank())
  
  return(plot)
}

# setup comparison along resource productivity gradient, define depth
setupComparison = function(mincoeff=0.1,maxcoeff=2,step=0.1,input_depth = 1500){
  coeff = seq(from=mincoeff, to=maxcoeff, by=step)
  stdzprod = 50
  y = list()
  
  groupvec = c("smallPel","largePel","demersals","totB")
  setuppal = c("#ADD8E6", "blue", "#8c8c8c", "#FFB6C1", "#FFD700")
  setupvec = c("setupBasicALT","setupBasic","setupBasic2",
               "setupVertical","setupVertical2")
  
  for (icoeff in 1:length(coeff)){
    dflux = 121+2.58*coeff[icoeff]*stdzprod
    basicbprod = 0.1*(dflux*(input_depth/150)^-0.86)
    
    if(basicbprod>=0.1*dflux) basicbprod=0.1*dflux
    p1=setupBasic(szprod = coeff[icoeff]*stdzprod, # small mesozoo production
                  lzprod = coeff[icoeff]*stdzprod, # large mesozoo production
                  bprodin = basicbprod, # benthos production
                  depth = input_depth) # Bottom depth
    p2=setupBasic2(szprod = coeff[icoeff]*stdzprod, # small mesozoo production
                   lzprod = coeff[icoeff]*stdzprod, # large mesozoo production
                   bprodin = basicbprod, # benthos production
                   depth = input_depth) # Bottom depth
    p3=setupVertical(szprod = coeff[icoeff]*stdzprod, # small mesozoo production
                     lzprod = coeff[icoeff]*stdzprod, # large mesozoo production
                     bprodin = basicbprod, # benthos production
                     depth = input_depth) # Bottom depth
    p4=setupVertical2(szprod = coeff[icoeff]*stdzprod,  # small mesozoo production
                      lzprod = coeff[icoeff]*stdzprod, # large mesozoo production
                      bprodin = basicbprod, # benthos production
                      depth = input_depth) # Bottom depth
    
    # alternative stable state
    p5=p1
    p5$u0[p5$ixFish]=1
    
    group_labels = c(p1$my_names[attr(p1$my_names, "names") %in%
                                   p1$groupnames[(p1$nResources + 1):length(p1$groupnames)]],
                     totB = "Total")
    
    for (i in 1:5){ # all setups
      assign(paste("sim", i, sep = ""),
             simulateFEISTY(p = get(paste("p", i, sep = "")),tEnd = 1000))
      Bpositive = get(paste("sim", i, sep = ""))$B
      Bpositive[Bpositive < 0] = 0
      p = get(paste("sim", i, sep = ""))$p
      sim = get(paste("sim", i, sep = ""))
      
      for (j in 1:p$nGroups) {
        y[[paste("sim", i, "coeff", coeff[icoeff], 
                 sep = "")]][[paste(p$groupnames[j + 4])]] = 
          sum(colMeans(Bpositive[round(0.6 * sim$nTime):sim$nTime,
                                 p$ix[[j]] - p$nResources]))
      }
      y[[paste("sim", i, "coeff", coeff[icoeff], sep = "")]][["totB"]] =
        sum(colMeans(Bpositive[round(0.6 * sim$nTime):sim$nTime,
                               p$ixFish - p$nResources]))
    }
  }
  
  filtered_lists1 <- y[grep("sim1.*", names(y))]
  filtered_lists2 <- y[grep("sim2.*", names(y))]
  filtered_lists3 <- y[grep("sim3.*", names(y))]
  filtered_lists4 <- y[grep("sim4.*", names(y))]
  filtered_lists5 <- y[grep("sim5.*", names(y))]
  
  # plots for biomass      
  ylim=c(0,max(unlist(y))+5)
  df=data.frame(x=coeff)
  ifig <- 1
  
  for (igroup in 1:4) {
    groupname=groupvec[igroup] 
    df[["setupBasicALT"]]  = unlist(sapply(filtered_lists5, 
                                           function(x)  x[[groupname]]))
    df[["setupBasic"]]     = unlist(sapply(filtered_lists1, 
                                           function(x) x[[groupname]]))
    df[["setupBasic2"]]    = unlist(sapply(filtered_lists2, 
                                           function(x) x[[groupname]]))
    df[["setupVertical"]]  = unlist(sapply(filtered_lists3, 
                                           function(x) x[[groupname]]))
    df[["setupVertical2"]] = unlist(sapply(filtered_lists4, 
                                           function(x) x[[groupname]]))
    if (igroup == 1){
      df[["setupVertical"]]=df[["setupVertical"]] + 
        unlist(sapply(filtered_lists3, function(x) x[["mesoPel"]]))
      df[["setupVertical2"]]=df[["setupVertical2"]] + 
        unlist(sapply(filtered_lists4, function(x) x[["mesoPel"]]))
    }
    if (igroup == 2){
      df[["setupVertical"]]=df[["setupVertical"]] + 
        unlist(sapply(filtered_lists3, function(x) x[["midwPred"]]))
      df[["setupVertical2"]]=df[["setupVertical2"]] + 
        unlist(sapply(filtered_lists4, function(x) x[["midwPred"]]))
    }
    
    df[["type"]]  = "biomass"
    df[['group']] = groupvec[igroup]
    
    df_long <- gather(df, key = "line", value = "y", -x,-type,-group)    
    if (ifig==1)   df_longlong=df_long
    if (!ifig==1)  df_longlong=rbind(df_longlong,df_long)
    plotorder=setupvec
    df_long$line <- factor(df_long$line, levels=plotorder)
    ifig=ifig+1
  }
  
  # plots for biomass ratio     
  ylim=ylim=c(0,1)
  df=data.frame(x=coeff)
  
  for (igroup in 1:3) { 
    groupname=groupvec[igroup] 
    df[["setupBasicALT"]]  = unlist(sapply(filtered_lists5, 
                                           function(x) x[[groupname]]/x[["totB"]]))
    df[["setupBasic"]]     = unlist(sapply(filtered_lists1, 
                                           function(x) x[[groupname]]/x[["totB"]]))
    df[["setupBasic2"]]    = unlist(sapply(filtered_lists2, 
                                           function(x) x[[groupname]]/x[["totB"]]))
    df[["setupVertical"]]  = unlist(sapply(filtered_lists3, 
                                           function(x) x[[groupname]]/x[["totB"]]))
    df[["setupVertical2"]] = unlist(sapply(filtered_lists4, 
                                           function(x) x[[groupname]]/x[["totB"]]))
    if (igroup == 1){
      df[["setupVertical"]]=df[["setupVertical"]] + 
        unlist(sapply(filtered_lists3, function(x) x[["mesoPel"]]/x[["totB"]]))
      df[["setupVertical2"]]=df[["setupVertical2"]] + 
        unlist(sapply(filtered_lists4, function(x) x[["mesoPel"]]/x[["totB"]]))
    }
    if (igroup == 2){
      df[["setupVertical"]]=df[["setupVertical"]] + 
        unlist(sapply(filtered_lists3, function(x) x[["midwPred"]]/x[["totB"]]))
      df[["setupVertical2"]]=df[["setupVertical2"]] + 
        unlist(sapply(filtered_lists4, function(x) x[["midwPred"]]/x[["totB"]]))
    }
    
    df[["type"]] = "ratio"
    df[['group']] = groupvec[igroup]
    
    df_long <- gather(df, key = "line", value = "y", -x,-type,-group)    
    df_longlong=rbind(df_longlong,df_long)
    plotorder=setupvec
    df_long$line <- factor(df_long$line, levels=plotorder)
    ifig=ifig+1
  }
  
  df_longlong$type <- factor(df_longlong$type, levels=c("biomass","ratio"))
  plotorder=groupvec
  df_longlong$group <- factor(df_longlong$group, levels=plotorder)
  plotorder=setupvec
  df_longlong$line <- factor(df_longlong$line, levels=plotorder)
  fig = ggplot(df_longlong, aes(x = x, y = y, color = line))+
    geom_line(linewidth = 0.7,alpha=0.9)+
    geom_point(size=1.5,alpha=0.9)+
    scale_color_manual(values = setuppal,
                       labels=c("setupBasic","setupBasicALT",
                                "setupBasic2","setupVertical","setupVertical2"),
                       breaks=c("setupBasic","setupBasicALT",
                                "setupBasic2","setupVertical","setupVertical2")) +
    labs(x = TeX("Zooplankton production (g m$^{-2}$ year$^{-1}$)"), 
         y = TeX("Relative to the total biomass                 Biomass (g m$^{-2}$)        "))+
    theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          axis.line = element_line(color = "black"),
          legend.title = element_blank(),
          legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.position = "bottom")+
    scale_x_continuous(breaks = c(0.1, 0.5, 1.0, 1.5, 2.0),
                       labels = c(0.1, 0.5, 1.0, 1.5, 2.0)*stdzprod)+
    facet_grid(rows=vars(type), cols=vars(group), scales = "free_y",
               labeller = labeller(group = as_labeller(group_labels))) +
    theme(strip.background = element_blank(), # edit theme for the facet labels
          strip.text = element_text(color = "blue"),
          strip.text.y = element_blank())
  
  grob <- ggplotGrob(fig)
  
  # Remove facets
  idx <- which(grob$layout$name %in% c("panel-2-4"))
  for (i in idx) grob$grobs[[i]] <- nullGrob()
  
  # Move x axes up
  # axis-b-4 needs to move up 2 rows
  idx <- which(grob$layout$name %in% c("axis-b-4"))
  grob$layout[idx, c("t", "b")] <- grob$layout[idx, c("t", "b")] - c(2);
  
  # Plot
  grid.newpage()
  grid.draw(grob)
  fig=grob
  return(fig)
}

# vary the number of stages 
analyseStagebiomass = function(nStages = c(3,6,9,12,15,18,21,24,27),vertical=F) {
  y = list()
  y2 = y
  ifig = 1
  
  # input preparation
  depths = c(100,1500) # shallow and deep
  zooprod = c(5,100) # low and high productivity
  vertdflux = 121 + 2.58 * zooprod # detrital flux out of photic zone low and high
  
  for(iprod in 1:length(zooprod)){ # low & high prod
    for (idepth in 1:length(depths)){ # shallow to deep water
      for (i in 1:length(nStages)) {
        
        if (vertical==F){
          basicbprod=0.1*(vertdflux[iprod]*(depths[idepth]/150)^-0.86)
          if(basicbprod>=0.1*vertdflux[iprod]) basicbprod=0.1*vertdflux[iprod]
          p=setupBasic2(szprod = zooprod[iprod], # small mesozoo production
                        lzprod = zooprod[iprod], # large mesozoo production
                        bprod = basicbprod,
                        depth  = depths[idepth],
                        nStages=nStages[i])
        }
        
        if (vertical==T){
          p=setupVertical2(szprod= zooprod[iprod], # small mesozoo production
                           lzprod = zooprod[iprod],# large mesozoo production
                           dfpho = vertdflux[iprod], # Detrital flux out photic zone
                           nStages=nStages[i],
                           depth=depths[idepth],
                           photic=150)
        }
        
        sim = simulateFEISTY(p=p)
        Bpositive=sim$B
        Bpositive[Bpositive<0]=0
        
        # all functional types
        for(j in 1:p$nGroups){
          y[[paste("Stage_", nStages[i], sep="")]][[paste(p$groupnames[j+4])]]= sum(colMeans(Bpositive[round(0.6*sim$nTime):sim$nTime,p$ix[[j]]-p$nResources])) 
        }
        
        # tot B     
        y[[paste("Stage_", nStages[i], sep="")]][["totB"]]= sum(colMeans(Bpositive[round(0.6*sim$nTime):sim$nTime,p$ixFish-p$nResources])) 
      }
      
      df=data.frame(x=nStages)
      for (i in 1:p$nGroups) {
        colname=p$groupnames[i+p$nResources]#paste0('y',i)
        df[[colname]]=unlist(sapply(y, function(x) x[[paste(p$groupnames[i+4])]]))
      }
      df[['totB']]=unlist(sapply(y, function(x) x[["totB"]]))
      df[['depth']]=depths[idepth]
      df[['zooprod']]=zooprod[iprod]
      
      df_long <- gather(df, key = "line", value = "y", -x,-zooprod,-depth)
      if (ifig==1)   df_longlong=df_long
      if (!ifig==1)  df_longlong=rbind(df_longlong,df_long)
      color_vec = c(p$my_palette[attr(p$my_palette, "names") %in% p$groupnames[(p$nResources+1):length(p$groupnames)]],"black")
      names(color_vec)[length(color_vec)] ="totB"
      name_vec = c(p$my_names[attr(p$my_names, "names") %in% p$groupnames[(p$nResources+1):length(p$groupnames)]],totB="Total")
      name_vec = unname(name_vec)
      
      ifig=ifig+1
      
    }
    
  }
  
  depth_labels=c(paste(depths[1],"m"),paste(depths[2],"m"))
  names(depth_labels)=depths
  zooprod_labels=c(paste(zooprod[1],"~","g ~ m^{-2} ~ year^{-1}"),paste(zooprod[2],"~","g ~ m^{-2} ~ year^{-1}"))
  names(zooprod_labels)=zooprod
  
  df_longlong$zooprod <- factor(df_longlong$zooprod, levels=c(zooprod[1],zooprod[2]))
  df_longlong$line <- factor(df_longlong$line, levels=c(p$groupnames[(p$nResources+1):length(p$groupnames)],"totB"))
  fig = ggplot(df_longlong, aes(x = x, y = y, color = line))+
    geom_line(linewidth = 0.7)+
    geom_point(size=1.5)+
    scale_color_manual(values = color_vec, labels=name_vec) +
    labs(x = "Stage number", y = TeX("Biomass (g m $^{-2}$)"))+
    theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          axis.line = element_line(color = "black"),
          legend.title = element_blank(), # remove legend title
          legend.key = element_rect(fill = "transparent", color = "transparent"), 
          legend.position = "bottom")+ 
    scale_x_continuous(breaks = unique(df_longlong$x))+
    facet_grid(rows=vars(zooprod), cols=vars(depth), scales = "free_y",
               labeller = labeller(zooprod = as_labeller(zooprod_labels,label_parsed),depth = depth_labels)) +
    theme(strip.background = element_blank(),
          strip.text = element_text(color = "blue"))
  
  return(fig)
}