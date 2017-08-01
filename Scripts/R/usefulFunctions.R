# # markers
#                   '<plot>',
#                   'type=text',
#                   'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/markers.txt',
#                   'color=black',
#                   'r1=0.99r',
#                   'r0=.80r',
#                   'show_links=yes',
#                   'link_dims=8p,8p,8p,4p,8p',
#                   'link_thickness=2p',
#                   'link_color=black',
#                   'label_size=13p',
#                   'label_font=default',
#                   'padding=2p',
#                   'rpadding=2p',
#                   'max_snuggle_distance=3r',
#                   'label_snuggle=yes',
#                   '</plot>',




# generador de huecolors
# 
hueGen<-function(n,from,to){
      x1<-round(seq(from,to,length.out=n),0)
      out1<-numeric()
      for (i in 1:length(x1)){
            out1[i]<-paste0('hue',ifelse(x1[i]<100,ifelse(x1[i]<10,paste0('00',x1[i]),paste0('0',x1[i])),x1[i]))
      }
      return(out1)
}