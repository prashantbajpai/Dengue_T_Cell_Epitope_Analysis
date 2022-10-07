library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(writexl)
library(patchwork)
library(stringr)
#source('/Users/prashant/R_themes/theme_publication.R')

read.and.filter.strain = function(seq_path){
  library(stringr)
  library(jsonlite)
  dat = fromJSON(seq_path)
  dat = dat[grep('polyprotein', dat$protName, ignore.case = T), ]
  dat = dat[!dat$date == 'N/A', ]
  dat = dat[!dat$date == "Inconsistent", ]
  dat$id_year = paste(dat$ncbiId, str_sub(dat$date, -4, -1), sep = ' ')
  dat$date[which(nchar(dat$date)==4)] = 
    paste('01', '01', dat$date[which(nchar(dat$date)==4)], sep = '/')
  dat$date[which(nchar(dat$date)==7)]=
    paste('01', dat$date[which(nchar(dat$date)==7)], sep = '/')
  print(names(table(nchar(dat$date))))
  dat$date = as.Date(dat$date, format = "%m/%d/%Y")
  dat = dat[dat$date > as.Date('2003', '%Y'), ]
  dat_list = split(dat, f = dat$virusType)
}

read.and.filter.epitopes = function(cd4_path, cd8_path, strain_list){
  epi_data_cd4 = read.csv(cd4_path)
  epi_data_cd4 = epi_data_cd4[c('Epitope.ID','aa_seq','start_new','end_new',
                                'Protein_name', 'Allele.Name', 'DV_status','E')]
  colnames(epi_data_cd4)[c(3,4)] = c('start','end')
  epi_data_cd8 = read.csv(cd8_path)
  epi_data_cd8 = epi_data_cd8[c('Epitope.ID','aa_seq','start_new','end_new',
                                'Protein_name', 'Allele.Name', 'DV_status','E')]
  colnames(epi_data_cd8)[c(3,4)] = c('start','end')
  epi_data_cd4_list = split(epi_data_cd4, f = epi_data_cd4$DV_status)
  epi_data_cd8_list = split(epi_data_cd8, f = epi_data_cd8$DV_status)
  
  epi_data_cd4_list = epi_data_cd4_list[names(strain_list)]
  epi_data_cd8_list = epi_data_cd8_list[names(strain_list)]
  return(list(epi_data_cd4_list = epi_data_cd4_list, 
              epi_data_cd8_list = epi_data_cd8_list))
}

align_matches = function(epi_data, sero_data, epi_id_col, sero_id_col,
                         epi_pepcol, sero_pepcol){
  library(Biostrings)
  mismatch_data = data.frame(matrix(data = NA, nrow = nrow(epi_data), 
                                    ncol = nrow(sero_data)))
  rownames(mismatch_data) = unlist(epi_data[epi_id_col])
  colnames(mismatch_data) = unlist(sero_data[sero_id_col])
  for(i in 1:nrow(epi_data)){
    for(j in seq_along(unlist(sero_data[sero_pepcol]))){
      matched_data = matchPattern(unlist(epi_data[epi_pepcol])[i], 
                                  unlist(sero_data[sero_pepcol])[j], 
                                  max.mismatch = 0, min.mismatch = 0)
      matched_num = nmismatch(unlist(epi_data[epi_pepcol])[i], matched_data)
      if(length(matched_num)==1){
        mismatch_data[i,j] = matched_num
      }
      if(length(matched_num)>1){
        mismatch_data[i,j] = sort(unique(matched_num))[1]
        print('wtf')
      }
    }
  }
  return(mismatch_data)
}

align_freq = function(df){
  library(reshape2)
  m0 = apply(df,1, function(x)sum(x==0,na.rm = T))
  test = data.frame(m0)
  test$m0 = (test$m0/ncol(df))*100
  test$step = cut(test$m0, breaks = c(-0.1,0,50,80,99,100), 
                  labels = c('0%','1%-50%','51%-80%','81%-99%','100%'))
  step_perc = table(test$step)
  step_perc = step_perc/nrow(df)*100
  step_perc = round(step_perc,2)
  step_perc = melt(step_perc)
  return(step_perc)
}

perc.plot = function(df, title){
  library(ggplot2)
  print(levels(df$Var1))
  mypalette = c('grey80', '#ffffb2', '#fecc5c', '#fd8d3c','#e31a1c')
  p1 = ggplot(df, aes(x = Serotype, y = value, fill = Var1)) + 
    geom_bar(stat = 'identity') + ylab('Percentage') + 
    xlab('Serotype') +
    #geom_text(aes(label = value), position = position_stack(vjust = 0.5),
    #          size = 2) +
    labs(title = title) + xlab('') +
    scale_y_continuous(labels = seq(0, 100, 10), breaks = seq(0, 100, 10)) +
    # scale_x_discrete(labels = c('DENV-1','DENV-2','DENV-3','DENV-4'),
    #                  drop = F) +
    scale_fill_manual(values = mypalette, name = '', 
                      labels = c('Epitopes not present \nin any isolates',
                                 'Epitopes present \nin 1%-50% of isolates',
                                 'Epitopes present \nin 51%-80% of isolates',
                                 'Epitopes present \nin 81%-99% of isolates',
                                 'Epitopes present \nin all isolates')) +
    theme_Publication() + 
    guides(fill=guide_legend(reverse = T, nrow = 1, byrow = T)) +
    #guides(fill=guide_legend(reverse = T)) +
    theme(legend.direction = 'vertical',
          axis.text.x = element_text(size = 16, angle = 90, vjust = 0.5),
          axis.text.y = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, face = 'bold', size = 20),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) 
  
  return(p1)
}

epi.metadata = function(df, epitope_data){
  library(reshape2)
  m0 = apply(df,1, function(x)sum(x==0,na.rm = T))
  test = data.frame(m0)
  test$m0 = (test$m0/ncol(df))*100
  test$step = cut(test$m0, breaks = c(-0.1,0,50,80,99,100), 
                  labels = c('0%','1%-50%','51%-80%','81%-99%','100%'))
  test$Epitope_ID = rownames(test)
  test$protein_name = rep('NA', nrow(test))
  test$aa_seq = rep('NA', nrow(test))
  test$start = rep('NA', nrow(test))
  test$end = rep('NA', nrow(test))
  test$Allele = rep('NA', nrow(test))
  test$Serotype = rep('NA', nrow(test))
  if(!grepl('start', paste(colnames(epitope_data), collapse = ','))){epitope_data$start = 'NA'}
  if(!grepl('end', paste(colnames(epitope_data), collapse = ','))){epitope_data$end = 'NA'}
  if(!grepl('DV_status', paste(colnames(epitope_data), collapse = ','))){epitope_data$DV_status = 'NA'}
  for(i in seq_along(epitope_data$Epitope.ID)){
    test$aa_seq[which(test$Epitope_ID %in% epitope_data$Epitope.ID[i])] =
      epitope_data$aa_seq[i]
    test$protein_name[which(test$Epitope_ID %in% epitope_data$Epitope.ID[i])] =
      epitope_data$Protein_name[i]
    test$start[which(test$Epitope_ID %in% epitope_data$Epitope.ID[i])] =
      epitope_data$start[i]
    test$end[which(test$Epitope_ID %in% epitope_data$Epitope.ID[i])] =
      epitope_data$end[i]
    test$Allele[which(test$Epitope_ID %in% epitope_data$Epitope.ID[i])] =
      epitope_data$Allele.Name[i]
    test$Serotype[which(test$Epitope_ID %in% epitope_data$Epitope.ID[i])] =
      epitope_data$DV_status[i]
  }
  test$protein_name[test$protein_name == 'Capsid'] = 'C'
  test$protein_name[test$protein_name == 'Env'] = 'E'
  return(test)
}

save.metadata = function(df, filname){
  library(xlsx)
  df = split(df, f = df$step)
  df = lapply(df, function(x){x = x[,-1]})
  filename = paste(filname, '_epitope_metadata.xlsx', sep = '')
  write.xlsx(df[[1]], filename, sheetName = names(df)[1], 
             row.names = F)
  write.xlsx(df[[2]], filename, sheetName = names(df)[2], 
             append = T, row.names = F)
  write.xlsx(df[[3]], filename, sheetName = names(df)[3], 
             append = T, row.names = F)
  write.xlsx(df[[4]], filename, sheetName = names(df)[4], 
             append = T, row.names = F)
  write.xlsx(df[[5]], filename, sheetName = names(df)[5], 
             append = T, row.names = F)
}

prot.dist.stepwise = function(df){
  df = df[c('step', 'protein_name')]
  df_list = split(df, f = df$step)
  prot_dist = lapply(df_list, function(x){
    plotdat = table(x$protein_name)
    plotdat = as.data.frame(plotdat)
    plotdat$perc = plotdat$Freq/sum(plotdat$Freq)*100
    plotdat$perc = round(plotdat$perc,2)
    return(plotdat)
  })
  return(prot_dist)
}

prot.dist.plot = function(df, backfill, title, x_order, x_label){
  df$Var1 = as.factor(df$Var1)
  x_order_fil = x_order[x_order %in% df$Var1]
  df$Var1 = factor(as.character(df$Var1), levels = x_order_fil)
  num_epi = sum(df$Freq)
  p1 = ggplot(df, aes(x = Var1, y = perc)) + 
    geom_bar(stat = 'identity', fill = 'black') + 
    #geom_text(aes(label = perc, y = perc+1)) +
    ylab('') +
    xlab('') + labs(title = title) +
    scale_x_discrete(limits = x_order, label = x_label) +
    theme(panel.background = element_rect(colour = NA, fill = backfill),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(colour="black"),
          axis.line.y = element_line(colour="black"),
          axis.ticks = element_line(),
          plot.background = element_rect(colour = NA, fill = NA),
          panel.border = element_rect(colour = NA, fill = NA),
          plot.title = element_text(face = 'bold', hjust = 0.5, size = 20),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 16, face = 'bold'),
          axis.text.y = element_text(size = 16)) +
    annotate('text', x = 9, y = max(df$perc, na.rm = T) + 5, 
             label = paste('Total: ', num_epi,sep = ''), size = 6)
  return(p1)
}

prot.dist.plot.compilation = function(x, title, prot_dist_data, yaxis = T){
  library(patchwork)
  library(ggpubr)
  mypalette = c('grey90', '#ffffb2', '#fecc5c', '#fd8d3c','#e83234')
  prot_list = prot_dist_data[[x]]
  p1_list = lapply(seq_along(prot_list), function(i){
    x_order = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a',
                'NS4b', 'NS5')
    x_label = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a',
                'NS4b', 'NS5')
    if(i==5){
      x = prot.dist.plot(prot_list[[i]], mypalette[i], 
                         names(prot_list)[i], x_order = x_order,
                         x_label = x_label)
    }else{
      x = prot.dist.plot(prot_list[[i]], mypalette[i], 
                         names(prot_list)[i], x_order = x_order,
                         x_label = NULL)
    }
  })
  if(yaxis){
    ytitle = text_grob('% of Epitopes', size = 18, rot = 90)
    p1 = {wrap_elements(ytitle) + wrap_plots(p1_list, nrow = 5) + 
        plot_layout(widths = c(0.03,1))} + 
      plot_annotation(title = title) & 
      theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5))
  }else{
    p1 = wrap_plots(p1_list, nrow = 5) + 
      plot_annotation(title = title) & 
      theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5))
  }
  return(p1)
}

color_table = function(df, E, Protein_name){
  library(dplyr)
  library(tidyr)
  library(janitor)
  library(formattable)
  colnames(df)[colnames(df) == E] = 'E'
  colnames(df)[colnames(df) == Protein_name] = 'Protein_name'
  df %>% 
    group_by(E) %>%
    count(Protein_name) %>%
    spread(Protein_name, n) %>%
    ungroup() %>%
    adorn_totals("col") %>%
    adorn_totals("row") -> test
  test = apply(test, 2,function(x){
    x[is.na(x)]=0
    return(x)
  })
  test = as.data.frame(test)
  test = test[c("E", "Capsid", "PreM", "Env", "NS1", "NS2a", "NS2b", "NS3",
                "NS4a", "NS4b", "NS5", "Total" )]
  # customblue = '#377eb8'
  # customblue0 = '#D7E5F0'
  # customGreen0 = "#DeF7E9"
  # customGreen = "#71CA97"
  # x = formattable(test, 
  #                 align =c("l","c","c","c","c", "c", "c", "c", "c", "c", 
  #                          "c", "c", "c"), 
  #                 list('E' = formatter("span", 
  #                                      style = ~ style(font.weight = "bold")),
  #                      area(1, col = c(-1,-12)) ~ color_tile(customblue0, 
  #                                                            customblue),
  #                      area(2, col = c(-1,-12)) ~ color_tile(customGreen0, 
  #                                                            customGreen)))
  return(test)
}

alignment_plots_data = function(x, epi_data, sero_data, metadata){
  if(nrow(x)>0){
    x$Epitope_ID = rownames(x)
    test = melt(x, id.vars = 'Epitope_ID')
    test$value = as.factor(as.character(test$value))
    test$date = rep(NA, nrow(test))
    class(test$date) = 'Date'
    test$id_year = rep('NA', nrow(test))
    test$protein_name = rep('NA', nrow(test))
    test$start = rep('NA', nrow(test))
    test$end = rep('NA', nrow(test))
    test$allele = rep('NA', nrow(test))
    test$aa_seq = rep('NA', nrow(test))
    test$step = rep('NA', nrow(test))
    test$virustype = rep('NA', nrow(test))
    test$dv_status = rep('NA', nrow(test))
    #test$megapool = rep('NA', nrow(test))
    if(!grepl('id_year', paste(colnames(sero_data), collapse = ','))){sero_data$id_year = 'NA'}
    
    for(k in 1: nrow(sero_data)){
      test$date[which(test$variable %in% sero_data$ncbiId[k])] =
        sero_data$date[k]
      test$id_year[which(test$variable %in% sero_data$ncbiId[k])] =
        sero_data$id_year[k]
      test$virustype[which(test$variable %in% sero_data$ncbiId[k])] =
        sero_data$virusType[k]
    }
    if(!grepl('start', paste(colnames(epi_data), collapse = ','))){epi_data$start = 'NA'}
    if(!grepl('end', paste(colnames(epi_data), collapse = ','))){epi_data$end = 'NA'}
    if(!grepl('DV_status', paste(colnames(epi_data), collapse = ','))){epi_data$DV_status = 'NA'}
    if(!grepl('step', paste(colnames(epi_data), collapse = ','))){epi_data$step = 'NA'}
    
    for(i in seq_along(epi_data$Epitope.ID)){
      test$protein_name[which(test$Epitope_ID %in% epi_data$Epitope.ID[i])] =
        epi_data$Protein_name[i]
      test$start[which(test$Epitope_ID %in% epi_data$Epitope.ID[i])] =
        epi_data$start[i]
      test$end[which(test$Epitope_ID %in% epi_data$Epitope.ID[i])] =
        epi_data$end[i]
      test$allele[which(test$Epitope_ID %in% epi_data$Epitope.ID[i])] =
        epi_data$Allele.Name[i]
      test$aa_seq[which(test$Epitope_ID %in% epi_data$Epitope.ID[i])] =
        epi_data$aa_seq[i]
      test$step[which(test$Epitope_ID %in% metadata$Epitope_ID[i])] =
        as.character(metadata$step[i])
      test$dv_status[which(test$Epitope_ID %in% epi_data$Epitope.ID[i])] =
        epi_data$DV_status[i]
    }
    if(!length(which(epi_data$DV_status == 'NA')) == nrow(epi_data)){
      test$dv_status = abbreviate(test$dv_status)
      test$dv_status = str_to_upper(test$dv_status)
      test$Epitope_ID = paste(test$Epitope_ID, ' (', test$protein_name, ' - ', 
                              test$dv_status,')', sep = '')
    }else{
      test$Epitope_ID = paste(test$Epitope_ID, ' (', test$protein_name, ')', sep = '')
    }
    plotdat = test[with(test, order(value, date)),]
  }else{
    plotdat = data.frame(matrix(ncol = 13, nrow = 0))
    colnames(plotdat) = c("Epitope_ID", "variable", "value", "date",
                          "id_year", "protein_name", "start", "end",
                          "allele", "aa_seq", "step", "virustype", "dv_status")
    plotdat = plotdat %>% mutate_all(as.character)
  }
  return(plotdat)
}

mytileplot = function(x, mytitle, xlim = c(-0.5, 5), ylim){
  x$virustype[is.na(x$value)] = NA
  x$step = factor(as.character(x$step), levels = c('0%', '1%-50%', '51%-80%', '81%-99%',
                                                   '100%'))
  x = x[order(x$step),]
  if(!grepl('id_year', paste(colnames(x), collapse = ','))){x$id_year = x$virustype}
  p1 = ggplot(x, aes(x=id_year,y=Epitope_ID)) + 
    geom_tile(aes(fill=virustype), width=.875, height=.875, color='white') +
    geom_label(aes(-0.3, Epitope_ID, label = Epitope_ID, fill = step), 
               size = 5, fontface = 'bold', show.legend = F) + 
    ylab('') + xlab('') + labs(title = mytitle) +
    scale_y_discrete(limits=rev(unique(x$Epitope_ID))) +
    theme_Publication() + 
    scale_fill_manual(name = '', 
                      values = c('0%' = 'grey80', '1%-50%' = '#ffffb2', 
                                 '51%-80%' = '#fecc5c', '81%-99%' = '#fd8d3c',
                                 '100%' = '#e04848', 
                                 "Japanese encephalitis virus" = '#377eb8',
                                 "Zika virus" = '#4daf4a',
                                 'SA14-JEV' = '#377eb8'),
                      breaks = levels(factor(x$virustype)),
                      na.value = 'grey90') +
    coord_cartesian(xlim = xlim, ylim = ylim, clip = 'off') +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, size = 16),
          plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = 16))
  return(p1)
}

epi.dist = function(df){
  df = split(df, f = df$virustype)
  perc = lapply(df, function(x){
    sero_count = nrow(x)
    x = table(x$step)
    x = as.data.frame(x)
    x$value = x$Freq/sum(x$Freq)*100
    x$sero_count = sero_count
    x = x[c('Var1', 'value', 'sero_count')]
    return(x)
  })
}

# perc.plot2 = function(df, title){
#   library(ggplot2)
#   mypalette = c('grey80', '#ffffb2', '#fecc5c', '#fd8d3c','#e31a1c')
#   p1 = ggplot(df, aes(x = Serotype, y = value, fill = Var1)) + 
#     geom_bar(stat = 'identity') + ylab('Percentage') + 
#     xlab('Serotype') +
#     #geom_text(aes(label = value), position = position_stack(vjust = 0.5),
#     #          size = 2) +
#     labs(title = title) +
#     scale_y_continuous(labels = seq(0, 100, 10), breaks = seq(0, 100, 10)) +
#     scale_x_discrete(drop = F) +
#     scale_fill_manual(values = mypalette, name = '', 
#                       labels = c('Epitopes not present in any isolates',
#                                  'Epitopes identical in 1%-50% of isolates',
#                                  'Epitopes identical in 51%-80% of isolates',
#                                  'Epitopes identical in 81%-99% of isolates',
#                                  'Epitopes identical in all isolates')) +
#     theme_Publication() + 
#     guides(fill=guide_legend(reverse = T, nrow = 3, byrow = T)) +
#     theme(legend.direction = 'vertical',
#           axis.text.x = element_text(size = 16, angle = 90, vjust = 0.5),
#           axis.text.y = element_text(size = 16),
#           plot.title = element_text(hjust = 0.5, face = 'bold', size = 20),
#           legend.key.size = unit(1, 'cm'),
#           legend.text = element_text(size = 16),
#           axis.title = element_text(size = 18)) 
  
#   return(p1)
# }

mytileplot_all_epi = function(x, color){
  y_order = x[c('Epitope_ID', 'protein_name')]
  y_order = y_order[order(y_order$protein_name),]
  x$Epitope_ID = factor(x$Epitope_ID, levels = unique(y_order$Epitope_ID))
  x$value = as.character(x$value)
  x$value[!is.na(x$value)] = x$serotype[!is.na(x$value)]
  x$value = factor(x$value, levels = c('DENV-1 Epitopes', 'DENV-2 Epitopes',
                                       'DENV-3 Epitopes', 'DENV-4 Epitopes'))
  label = unique(x$protein_name)
  if(label %in% c('NS3', 'NS5')){
    p1 = ggplot(x, aes(x=variable,y=Epitope_ID, fill = value)) + 
      geom_tile(aes(fill=value), width=.875, height=.875, color='#f0f0f0') +
      scale_x_discrete(limits=unique(x$variable)) +
      #scale_y_discrete(limits=unique(y_order$Epitope_ID)) +
      theme_Publication() + 
      ylab('Epitope ID') + xlab('') +
      scale_fill_manual(values = color, name = '',
                        breaks = levels(x$value),
                        na.value = 'grey90') +
      theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 1),
            axis.text.y = element_text(size = 8),
            legend.key = element_rect(size = 2, color = 'white'),
            legend.key.size = unit(1.5, 'lines'),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            legend.position = 'none') 
  }else{
    p1 = ggplot(x, aes(x=variable,y=Epitope_ID, fill = value)) + 
      geom_tile(aes(fill=value), width=.875, height=.875, color='#f0f0f0') +
      scale_x_discrete(limits=unique(x$variable)) +
      #scale_y_discrete(limits=unique(y_order$Epitope_ID)) +
      theme_Publication() + 
      ylab('Epitope ID') + xlab('') +
      scale_fill_manual(values = color, name = '',
                        breaks = levels(x$value),
                        na.value = 'grey90') +
      theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 1),
            axis.text.y = element_text(size = 8),
            legend.key = element_rect(size = 2, color = 'white'),
            legend.key.size = unit(1.5, 'lines'),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            legend.position = 'none') 
  }
}

alignment_plots_data2 = function(x, epi_data, sero_data, metadata){
  if(nrow(x)>0){
    x$Epitope_ID = rownames(x)
    test = melt(x, id.vars = 'Epitope_ID')
    test$value = as.factor(as.character(test$value))
    test$date = rep(NA, nrow(test))
    class(test$date) = 'Date'
    test$id_year = rep('NA', nrow(test))
    test$protein_name = rep('NA', nrow(test))
    test$start = rep('NA', nrow(test))
    test$end = rep('NA', nrow(test))
    test$allele = rep('NA', nrow(test))
    test$aa_seq = rep('NA', nrow(test))
    test$step = rep('NA', nrow(test))
    test$virustype = rep('NA', nrow(test))
    test$dv_status = rep('NA', nrow(test))
    #test$megapool = rep('NA', nrow(test))
    if(!grepl('id_year', paste(colnames(sero_data), collapse = ','))){sero_data$id_year = 'NA'}
    
    for(k in 1: nrow(sero_data)){
      test$date[which(test$variable %in% sero_data$ncbiId[k])] =
        sero_data$date[k]
      test$id_year[which(test$variable %in% sero_data$ncbiId[k])] =
        sero_data$id_year[k]
      test$virustype[which(test$variable %in% sero_data$ncbiId[k])] =
        sero_data$virusType[k]
    }
    if(!grepl('start', paste(colnames(epi_data), collapse = ','))){epi_data$start = 'NA'}
    if(!grepl('end', paste(colnames(epi_data), collapse = ','))){epi_data$end = 'NA'}
    if(!grepl('DV_status', paste(colnames(epi_data), collapse = ','))){epi_data$DV_status = 'NA'}
    if(!grepl('step', paste(colnames(epi_data), collapse = ','))){epi_data$step = 'NA'}
    
    for(i in seq_along(epi_data$Epitope.ID)){
      test$protein_name[which(test$Epitope_ID %in% epi_data$Epitope.ID[i])] =
        epi_data$Protein_name[i]
      test$start[which(test$Epitope_ID %in% epi_data$Epitope.ID[i])] =
        epi_data$start[i]
      test$end[which(test$Epitope_ID %in% epi_data$Epitope.ID[i])] =
        epi_data$end[i]
      test$allele[which(test$Epitope_ID %in% epi_data$Epitope.ID[i])] =
        epi_data$Allele.Name[i]
      test$aa_seq[which(test$Epitope_ID %in% epi_data$Epitope.ID[i])] =
        epi_data$aa_seq[i]
      test$step[which(test$Epitope_ID %in% metadata$Epitope_ID[i])] =
        as.character(metadata$step[i])
      test$dv_status[which(test$Epitope_ID %in% epi_data$Epitope.ID[i])] =
        epi_data$DV_status[i]
    }
    if(!length(which(epi_data$DV_status == 'NA')) == nrow(epi_data)){
      test$dv_status = abbreviate(test$dv_status)
      test$dv_status = str_to_upper(test$dv_status)
      test$Epitope_ID = paste(test$Epitope_ID, '(', test$dv_status,')', sep = '')
      #test$variable = paste(test$variable, '\n (', test$dv_status,')', sep = '')
    }else{
      # test$Epitope_ID = paste(test$Epitope_ID, ' (', test$protein_name, ')', sep = '')
    }
    plotdat = test[with(test, order(dv_status, date)),]
  }else{
    plotdat = data.frame(matrix(ncol = 13, nrow = 0))
    colnames(plotdat) = c("Epitope_ID", "variable", "value", "date",
                          "id_year", "protein_name", "start", "end",
                          "allele", "aa_seq", "step", "virustype", "dv_status")
    plotdat = plotdat %>% mutate_all(as.character)
  }
  return(plotdat)
}

barcode_like_plot = function(meta){
  meta = split(meta, f = meta$Serotype)
  
  #print(names(meta))
  serowise = lapply(meta, function(serodata){
    serodata$protein_name = factor(serodata$protein_name, 
                                   levels = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b',
                                              'NS3', 'NS4a', 'NS4b', 'NS5'))
    serodata = serodata[order(serodata$protein_name),]
    serodata$xlab = seq_len(nrow(serodata))
    serodata$xlab = factor(serodata$xlab, levels = serodata$xlab)
    steps = split(serodata, f = serodata$step)
    
    rectdat = split(serodata, f = serodata$protein_name)
    rectdat = lapply(rectdat, function(x){
      if(nrow(x)>0){
        x$xlab = as.numeric(as.character(x$xlab))
        out = data.frame(xmin = min(x$xlab), xmax = max(x$xlab) + 1)
      }
    })
    rectdat = bind_rows(rectdat, .id = 'protein')
    rectdat$xmin[1] = 0
    rectdat$ymin = -Inf
    rectdat$ymax = Inf
    rectdat$protein = factor(rectdat$protein, 
                             levels = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b',
                                        'NS3', 'NS4a', 'NS4b', 'NS5'))
    barwidth = (max(rectdat$xmax)/395)*0.8
    print(names(steps))
    step_wise = lapply(steps, function(x){
      x$ylab = 1
      bordercol = unique(as.character(x$step))
      bordercol = ifelse(bordercol == '0%', 'grey50', 
                         ifelse(bordercol == '1%-50%', '#ffffb2',
                                ifelse(bordercol == '51%-80%', '#fecc5c',
                                       ifelse(bordercol == '81%-99%', '#fd8d3c',
                                              '#e31a1c'))))
      bound_dat = data.frame(xmin = min(rectdat$xmin), 
                             xmax = max(rectdat$xmax),
                             ymin=-Inf,ymax=Inf)
      p_out = ggplot(x, aes(x = xlab, y = ylab, fill = protein_name)) +
        
        # geom_rect(data = bound_dat, inherit.aes = F,
        #           aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
        #           alpha = 1,
        #           fill = 'red',
        #           size = 3) +
        geom_rect(data = rectdat, inherit.aes = F, 
                  aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
                      fill = protein), alpha = 0.2,
                  size = 1.5, color = NA) +
        geom_bar(stat = 'identity', width = barwidth) +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_discrete(drop = F) +
        scale_fill_manual(values = c("#6BAED6", "#fb9a99",
                                     "#525252", '#6A3D9A',
                                     "#FD8D3C", "#FDD0A2",
                                     "#02818a", '#ae017e', '#f768a1', "#8c510a"),
                          name = '', drop = F) +
        theme_void() +
        theme(legend.position = 'none',
              plot.margin = unit(c(0, 0, 0, 0), "lines"),
              plot.background = element_rect(fill = 'white', colour = NA))
    })
    # p_final = step_wise[[1]] + step_wise[[2]] + step_wise[[3]] +
    #   step_wise[[4]] + step_wise[[5]] + plot_layout(ncol = 1, guides = 'collect') &
    #   theme(legend.position = 'none')
    
  })
}

perc.plot2 = function(df, title = '', legendpos = 'none'){
  library(ggplot2)
  mypalette = c('grey80', '#ffffb2', '#fecc5c', '#fd8d3c','#e31a1c')
  #mypalette = tail(mypalette, length(unique(df$Var1)))
  df$Var1 = factor(as.character(df$Var1), 
                   levels = c("0%", "1%-50%", "51%-80%", 
                              "81%-99%", "100%"))
  p1 = ggplot(df, aes(x = Serotype, y = value, fill = Var1)) + 
    geom_bar(stat = 'identity') + ylab('Percentage') + 
    xlab('Serotype') +
    #geom_text(aes(label = value), position = position_stack(vjust = 0.5),
    #          size = 2) +
    labs(title = title) + 
    xlab('') +
    scale_y_continuous(labels = seq(0, 100, 10), breaks = seq(0, 100, 10), limits = c(0, 100.5)) +
    scale_x_discrete(drop = F) +
    scale_fill_manual(values = mypalette, name = '', 
                      labels = c('Epitopes not present in any isolates',
                                 'Epitopes present in 1%-50% of isolates',
                                 'Epitopes present in 51%-80% of isolates',
                                 'Epitopes present in 81%-99% of isolates',
                                 'Epitopes present in all isolates'), drop = F) +
    theme_Publication() + 
    guides(fill=guide_legend(reverse = T, nrow = 3, byrow = T)) +
    #guides(fill=guide_legend(reverse = T)) +
    theme(legend.direction = 'vertical',
          axis.text.x = element_text(size = 16, angle = 90, vjust = 0.5, face = 'bold'),
          axis.text.y = element_text(size = 16, face = 'bold'),
          plot.title = element_text(hjust = 0.5, face = 'bold', size = 20),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          legend.position = legendpos,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) 
  
  return(p1)
}

bar_like_plot = function(meta, legendpos = 'none'){
  meta = split(meta, f = meta$Serotype)
  serowise = lapply(meta, function(serodata){
    serodata$protein_name = factor(serodata$protein_name, 
                                   levels = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b',
                                              'NS3', 'NS4a', 'NS4b', 'NS5'))
    protwise = split(serodata, f = serodata$protein_name)
    protwise = lapply(seq_along(protwise), function(i){
      protdat = protwise[[i]]
      protdat = data.frame(table(protdat$step))
      protdat$protein_name = names(protwise)[i]
      protdat$perc = protdat$Freq/sum(protdat$Freq) * 100
      sumdat = data.frame('sum' = sum(protdat$Freq), 
                          'protein_name' = names(protwise)[i])
      return(list(protdat, sumdat))
    })
    protwise_perc = lapply(protwise, function(x){x[[1]]})
    protwise_perc = bind_rows(protwise_perc)
    protwise_total = lapply(protwise, function(x){x[[2]]})
    protwise_total = bind_rows(protwise_total)
    protwise_perc$protein_name = factor(protwise_perc$protein_name, 
                                        levels = c('C', 'PreM', 'E', 'NS1', 
                                                   'NS2a', 'NS2b', 'NS3', 
                                                   'NS4a', 'NS4b', 'NS5'))
    protwise_perc$Var1 = factor(as.character(protwise_perc$Var1), 
                                levels = c("0%", "1%-50%", "51%-80%", 
                                           "81%-99%", "100%"))
    p1 = ggplot(protwise_perc, aes(x = protein_name, y = perc, fill = Var1)) +
      geom_bar(stat = 'identity') +
      geom_text(data = protwise_total, aes(x = protein_name, y = 103, label = sum),
                inherit.aes = F, fontface = 'bold', size = 6) +
      theme_Publication() +
      scale_fill_manual(values = c('grey80', '#ffffb2', '#fecc5c', '#fd8d3c',
                                   '#e31a1c')) +
      scale_y_continuous(limits = c(0, 103), labels = seq(0, 100, 10), 
                         breaks = seq(0, 100, 10)) +
      xlab('') + ylab('Percentage') +
      theme(legend.position = legendpos,
            axis.text.x = element_text(face = 'bold', colour = 'black', size = 14),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 18, face = 'bold'),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  })
}

bar_like_plot2 = function(meta, legendpos = 'none'){
  meta = split(meta, f = meta$virustype)
  serowise = lapply(meta, function(serodata){
    serodata$protein_name = gsub('Capsid', 'C', serodata$protein_name)
    serodata$protein_name = gsub('Env', 'E', serodata$protein_name)
    serodata$protein_name = factor(serodata$protein_name, 
                                   levels = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b',
                                              'NS3', 'NS4a', 'NS4b', 'NS5'))
    protwise = split(serodata, f = serodata$protein_name)
    protwise = lapply(seq_along(protwise), function(i){
      protdat = protwise[[i]]
      if(nrow(protdat) > 0){
        protdat = data.frame(table(protdat$step))
        protdat$protein_name = names(protwise)[i]
        protdat$perc = protdat$Freq/sum(protdat$Freq) * 100
        sumdat = data.frame('sum' = sum(protdat$Freq), 
                            'protein_name' = names(protwise)[i])
        return(list(protdat, sumdat))
      }
    })
    protwise_perc = lapply(protwise, function(x){x[[1]]})
    protwise_perc = bind_rows(protwise_perc)
    protwise_total = lapply(protwise, function(x){x[[2]]})
    protwise_total = bind_rows(protwise_total)
    protwise_perc$protein_name = factor(protwise_perc$protein_name, 
                                        levels = c('C', 'PreM', 'E', 'NS1', 
                                                   'NS2a', 'NS2b', 'NS3', 
                                                   'NS4a', 'NS4b', 'NS5'))
    protwise_perc$Var1 = factor(as.character(protwise_perc$Var1), 
                                levels = c("0%", "1%-50%", "51%-80%", 
                                           "81%-99%", "100%"))
    p1 = ggplot(protwise_perc, aes(x = protein_name, y = perc, fill = Var1)) +
      geom_bar(stat = 'identity') +
      geom_text(data = protwise_total, aes(x = protein_name, y = 105, label = sum),
                inherit.aes = F, fontface = 'bold', size = 6) +
      theme_Publication() +
      scale_fill_manual(values = c('grey80', '#ffffb2', '#fecc5c', '#fd8d3c',
                                   '#e31a1c'), drop = F) +
      scale_y_continuous(limits = c(0, 105), labels = seq(0, 100, 10), 
                         breaks = seq(0, 100, 10)) +
      scale_x_discrete(drop = F) +
      xlab('') + ylab('Percentage') +
      theme(legend.position = legendpos,
            axis.text.x = element_text(face = 'bold', colour = 'black', size = 14),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 18, face = 'bold'),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  })
}

barcode_labels = function(meta){
  meta = split(meta, f = meta$Serotype)
  
  #print(names(meta))
  serowise = lapply(meta, function(serodata){
    serodata$protein_name = factor(serodata$protein_name, 
                                   levels = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b',
                                              'NS3', 'NS4a', 'NS4b', 'NS5'))
    serodata = serodata[order(serodata$protein_name),]
    serodata$xlab = seq_len(nrow(serodata))
    serodata$xlab = factor(serodata$xlab, levels = serodata$xlab)
    steps = split(serodata, f = serodata$step)
    
    rectdat = split(serodata, f = serodata$protein_name)
    rectdat = lapply(rectdat, function(x){
      if(nrow(x)>0){
        x$xlab = as.numeric(as.character(x$xlab))
        out = data.frame(xmin = min(x$xlab), xmax = max(x$xlab) + 1)
      }
    })
    rectdat = bind_rows(rectdat, .id = 'protein')
    rectdat$xmin[1] = 0
    rectdat$ymin = -Inf
    rectdat$ymax = Inf
    rectdat$protein = factor(rectdat$protein, 
                             levels = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b',
                                        'NS3', 'NS4a', 'NS4b', 'NS5'))
    barwidth = (max(rectdat$xmax)/395)*0.8
    print(names(steps))
    step_wise = lapply(steps, function(x){
      x$ylab = 1
      bordercol = unique(as.character(x$step))
      bordercol = ifelse(bordercol == '0%', 'grey50', 
                         ifelse(bordercol == '1%-50%', '#ffffb2',
                                ifelse(bordercol == '51%-80%', '#fecc5c',
                                       ifelse(bordercol == '81%-99%', '#fd8d3c',
                                              '#e31a1c'))))
      bound_dat = data.frame(xmin = min(rectdat$xmin), 
                             xmax = max(rectdat$xmax),
                             ymin=-Inf,ymax=Inf)
      p_out = ggplot(x, aes(x = xlab, y = ylab, fill = protein_name)) +
        
        # geom_rect(data = bound_dat, inherit.aes = F,
        #           aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
        #           alpha = 1,
        #           fill = 'red',
        #           size = 3) +
        geom_rect(data = rectdat, inherit.aes = F, 
                  aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
                      fill = protein), alpha = 0.8,
                  size = 1.5, color = NA) +
        #geom_bar(stat = 'identity', width = barwidth) +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_discrete(drop = F) +
        scale_fill_manual(values = c("#6BAED6", "#fb9a99",
                                     "#525252", '#6A3D9A',
                                     "#FD8D3C", "#FDD0A2",
                                     "#02818a", '#ae017e', '#f768a1', "#8c510a"),
                          name = '', drop = F) +
        theme_void() +
        theme(legend.position = 'none',
              plot.margin = unit(c(0, 0, 0, 0), "lines"),
              plot.background = element_rect(fill = 'white', colour = NA))
    })
    # p_final = step_wise[[1]] + step_wise[[2]] + step_wise[[3]] +
    #   step_wise[[4]] + step_wise[[5]] + plot_layout(ncol = 1, guides = 'collect') &
    #   theme(legend.position = 'none')
    
  })
}

vaccine_post_alignment_prcossing = function(align_dat, vaccine_name, vaccine_df, epi_data, epi_metadata){
  align_dat = lapply(align_dat, function(df){
    df$rowsum = apply(df,1, function(x)sum(x==0,na.rm = T))
    df = df[df$rowsum!=0,]
    df$rowsum = NULL
    return(df)
  })

  align_plotdat = lapply(align_dat, function(df){
    dat = alignment_plots_data(df, epi_data,
                               vaccine_df, epi_metadata)
  })
  align_plotdat = lapply(align_plotdat, function(x){
    x$date = as.character(x$date)
    return(x)
  })
  align_plotdat = bind_rows(align_plotdat)

  align_plotdat = align_plotdat[c('Epitope_ID', 'step','protein_name',
                                                'virustype', 'value', 'dv_status')]
  align_plotdat$dv_status = vaccine_name
  align_plotdat = align_plotdat[!is.na(align_plotdat$value),]

  align_dist = epi.dist(align_plotdat)
  align_dist = bind_rows(align_dist, .id = 'Serotype')
  align_dist = align_dist[!align_dist$Var1 == 'NA',]
  align_dist$Serotype = paste(align_dist$Serotype, ' (', 
                                         align_dist$sero_count, ')', 
                                         sep = '')
  xlabels = unique(align_dist$Serotype)
  align_dist$Serotype = factor(align_dist$Serotype, 
                                          levels = xlabels)
  return(align_dist)
}

vaccine_metadata = function(align_dat, vaccine_name, vaccine_df, epi_data, epi_metadata){
  align_dat = lapply(align_dat, function(df){
    df$rowsum = apply(df,1, function(x)sum(x==0,na.rm = T))
    df = df[df$rowsum!=0,]
    df$rowsum = NULL
    return(df)
  })
  
  align_plotdat = lapply(align_dat, function(df){
    dat = alignment_plots_data(df, epi_data,
                               vaccine_df, epi_metadata)
  })
  align_plotdat = lapply(align_plotdat, function(x){
    x$date = as.character(x$date)
    return(x)
  })
  align_plotdat = bind_rows(align_plotdat)
  align_plotdat$dv_status = vaccine_name
  align_plotdat = align_plotdat[!is.na(align_plotdat$value),]
  align_plotdat$Epitope_ID = gsub(' .*', '', align_plotdat$Epitope_ID)
  return(align_plotdat)
}

countrywise_alignment = function(contry_list, cd4_list, cd8_list){
  cd4_list = cd4_list[names(contry_list)]
  cd8_list = cd8_list[names(contry_list)]
  print('Perfoming alignment of CD4 epitopes')
  cd4_alignment = lapply(seq_along(cd4_list), function(i){
    match_data = align_matches(epi_data = cd4_list[[i]],
                               sero_data = contry_list[[i]],
                               epi_id_col = 'Epitope.ID',
                               sero_id_col = 'ncbiId',
                               epi_pepcol = 'aa_seq',
                               sero_pepcol = 'aaSeq')
    return(match_data)
  })
  names(cd4_alignment) = names(contry_list)
  print('Perfoming alignment of CD8 epitopes')
  cd8_alignment = lapply(seq_along(cd8_list), function(i){
    match_data = align_matches(epi_data = cd8_list[[i]],
                               sero_data = contry_list[[i]],
                               epi_id_col = 'Epitope.ID',
                               sero_id_col = 'ncbiId',
                               epi_pepcol = 'aa_seq',
                               sero_pepcol = 'aaSeq')
    return(match_data)
  })
  names(cd8_alignment) = names(contry_list)
  print('Calculating aligment frequencies and writing plots')
  cd4_perc = lapply(cd4_alignment, align_freq)
  cd4_perc = bind_rows(cd4_perc, .id = 'Serotype')
  cd8_perc = lapply(cd8_alignment, align_freq)
  cd8_perc = bind_rows(cd8_perc, .id = 'Serotype')
  # cd4_perc$Serotype = factor(cd4_perc$Serotype, levels = c('Dengue virus 1', 
  #                                                          'Dengue virus 2',
  #                                                          'Dengue virus 3',
  #                                                          'Dengue virus 4'))
  # cd8_perc$Serotype = factor(cd8_perc$Serotype, levels = c('Dengue virus 1', 
  #                                                          'Dengue virus 2',
  #                                                          'Dengue virus 3',
  #                                                          'Dengue virus 4'))
  # p1 = perc.plot2(cd4_perc, 'CD4 Epitopes')
  # p2 = perc.plot2(cd8_perc, 'CD8 Epitopes')
  # p_final = p1 + p2 + plot_layout(ncol = 2, guides = 'collect') &
  #   theme(legend.position = 'none')
  # return(p_final)
  out = list('CD4_perc' = cd4_perc, 'CD8_perc' = cd8_perc, 'CD4_align' = cd4_alignment, 'cd8_alignment' = cd8_alignment)
}

custom_upset = function(metadf, minsetsize, left_ylim = 300){
  library(ComplexHeatmap)
  library(ComplexUpset)
  library(ggbeeswarm)
  library(ggpattern)

  upsetdat = metadf[c('Epitope_ID', 'Country', 'protein_name')]
  listset = split(upsetdat, f = upsetdat$Country)
  listset = lapply(listset, function(x){x = x$Epitope_ID})
  upsetmat = list_to_matrix(listset)
  upsetmat = as.data.frame(upsetmat)
  upsetmat$Epitope_ID = rownames(upsetmat)

  upsetdat = upsetdat[!duplicated(upsetdat$Epitope_ID),]
  upsetdat = upsetdat[c('Epitope_ID', 'protein_name')]
  plotdat = merge(upsetmat, upsetdat, by = 'Epitope_ID')
  plotdat$protein_name[plotdat$protein_name == 'E'] = 'Env'
  plotdat$protein_name[plotdat$protein_name == 'C'] = 'Capsid'
  plotdat$protein_name = factor(plotdat$protein_name, 
                                levels = c('Capsid', 'PreM', 'Env', 'NS1', 'NS2a',
                                           'NS2b', 'NS3', 'NS4a', 'NS4b', 'NS5'))

  #colnames(plotdat) = str_to_upper(colnames(plotdat))
  #colnames(plotdat)[2:7] = c("DENVAX", "DPIV", "TV003", "DENGVAXIA", "TDEN", "TVDV")

  combination = setdiff(colnames(plotdat), c('Epitope_ID', 'protein_name'))
  print(combination)
  #set.seed(100)
  p1 = upset(plotdat, combination, width_ratio = 0.2,
             name = 'Country',
             themes = upset_modify_themes(
               list(
                 'intersections_matrix'=theme(text=element_text(size=18, face = 'bold'),
                                              axis.text = element_text(colour = 'black'))
               )),
             # matrix = (
             #   intersection_matrix(
             #     geom = geom_point()
             #   )
             #   + scale_y_discrete(position = 'right')
             # ),
             min_size = minsetsize,
             sort_intersections='descending',
             base_annotations = list(
               'Intersection size' = intersection_size(
                 text = list(size = 6, color = 'black', fontface = 'bold'),
                 mapping = aes(fill='bars_fill', color = 'bar_color'))
               + ylab('Epitopes present \nin all strains')
               + scale_fill_manual(values = c('bars_fill'='grey70'), guide = 'none')
               #+ scale_color_manual(values = c('bar_color' = 'black'), guide = 'none')
               + theme(axis.title.y = element_text(face = 'bold',
                                                   size = 18),
                       axis.text.y = element_text(size = 16, 
                                                  colour = 'black', face = 'bold'))
             ),
             set_sizes = upset_set_size(aes(width = 0.1)) + 
               ylab('Epitopes present \nin all strains') 
             + geom_text(aes(label=..count..), hjust=1.1, stat='count', 
                         fontface = 'bold', size = 6)
             #+ geom_bar_pattern(aes(pattern = ))
             + expand_limits(y=left_ylim)
             + theme(axis.text.x = element_blank(), 
                     axis.title = element_text(size = 18, face = 'bold'),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()),
             sort_sets = 'descending',
             # queries = list(
             #   upset_query(set = 'DENVAX', fill = '#e4c97d'),
             #   upset_query(set = 'DPIV', fill = '#fb6a4a'),
             #   upset_query(set = 'NIH', fill = '#fb6a4a'),
             #   upset_query(set = 'SEN', fill = '#e4c97d'),
             #   upset_query(set = 'TDEN', fill = '#fb6a4a'),
             #   upset_query(set = 'TVDV', fill = '#e4c97d')),
             annotations = list(
               'Protein Distribution'=(
                 ggplot(mapping = aes(fill = protein_name))
                 + geom_bar(stat='count', position=position_fill(reverse = T), 
                            na.rm=TRUE)
                 # + geom_text(aes(label=!!aes_percentage(relative_to='intersection')),
                 #             stat='count',
                 #             position=position_fill(vjust = 0.5, reverse = T), size = 5,
                 #             fontface = 'bold')
                 + ylab('Protein Distribution')
                 + scale_y_continuous(labels=scales::percent_format())
                 + scale_fill_manual(values = c("#6BAED6", "#fb9a99",
                                                "#525252", '#6A3D9A',
                                                "#FD8D3C", "#FDD0A2",
                                                "#02818a", '#ae017e', 
                                                '#f768a1', "#8c510a"), 
                                     name = 'Dengue \npolyprotein', drop = F)
                 #+ scale_color_manual(values=c('show'='black', 'hide'='transparent'), 
                 #                    guide=FALSE)
                 + theme(axis.title.y = element_text(face = 'bold', size = 18),
                         axis.text.y = element_text(size = 16, colour = 'black',
                                                    face = 'bold'),
                         legend.text = element_text(size = 16),
                         legend.key.size = unit(1, 'cm'),
                         legend.title = element_text(size = 16),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank())
               )
             )
  )

}


