options(stringsAsFactors = F)
if (interactive()) {
  setwd('~/d/sci/src/prp_binder_screening')
}
library(sqldf)

# imgmode = 'png'
imgmode = 'png'
imgsave = get(imgmode) # get the function that saves as that type of image
if (imgmode == 'png') {
  resx = 600 # multiplier for image width, height, and resolution
} else if (imgmode == 'pdf') {
  resx = 1
}

percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='f'),"%",sep="") ) )
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}


# parameters for table of libraries screened
liborder = data.frame(library=c('f19','chiral','cdot1','chembridge','cdot2'),
                      lorder=c(1,4,2,5,3),
                      disp=c('Broad 19F library','Schreiber chiral collection','Broad 1st gen STD library','Chembridge high-solubility set','Broad 2nd gen STD library'),
                      color=rev(c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e')))

  
frags = read.table('data/fragments.tsv', sep='\t', header=T, quote='', comment.char='')

frags$color = liborder$color[match(frags$library, liborder$library)]

# check descriptors - make sure all are present?
sum(is.na(frags$hba))
sum(is.na(frags$hbd))
sum(is.na(frags$mw))
sum(is.na(frags$logp))
sum(is.na(frags$rotb))
sum(is.na(frags$tpsa))
range(frags$logp)

#### TABLE 1

table1 = sqldf("
               select   l.disp, count(*) n, sum(hit) hits, sum(std_attempted) stds, sum(trosy_attempted) trosys, sum(trosy_validated) validated
               from     frags f, liborder l
               where    f.library = l.library
               group by 1
               order by l.lorder
               ;")

sum(table1$n)

table1$color = rev(c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e'))

table1$wrap = gsub('(.{1,20})(\\s|$)', '\\1\n', table1$disp) # https://stackoverflow.com/a/2352006/3806692
table1$wrap = trimws(table1$wrap) # remove terminal line breaks

resx=300
#### FIGURE 1 - properties of library
imgsave(paste('display_items/figure-1.',imgmode,sep=''),width=6.5*resx,height=2.75*resx,res=resx)

layout_matrix = matrix(c(1,2,3,1,2,4),nrow=2,byrow=T)
layout(layout_matrix)

par(mar=c(4,4,3,1))

# stacked bar
# barplot(as.matrix(table1[,c('n')]), col=table1$color, beside=FALSE, border=NA, axes=F)
table1$cumn = cumsum(table1$n)
table1$mid = table1$cumn - 0.5*table1$n
plot(NA, NA, xlim=c(0,1.5), ylim=c(0,7000), xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=c(0,1.5), labels=NA, lwd=1, lwd.ticks=0)
axis(side=2, at=(0:7)*1000, labels=formatC((0:7)*1000,big.mark=','), lwd=1, lwd.ticks=1, las=2)
for (i in nrow(table1):1) {
  rect(xleft=0.05, xright=0.5, ybottom=table1$cumn[i] - table1$n[i], ytop=table1$cumn[i], col=table1$color[i], border=NA)
  text(x=0.5, y=table1$mid[i], pos=4, labels=table1$wrap[i], col=table1$color[i], font=2, cex=0.9)
}
mtext('A', side=3, cex=2, adj = -0.3, line = 0.5)

compliant_color = '#FF9912'

# Congreve 2003 - Rule of 3 - https://www.ncbi.nlm.nih.gov/pubmed/14554012
# overall Ro3 compliance:
mean(frags$mw <= 300 & frags$logp <= 3 & frags$hba <= 3 & frags$hbd <= 3)

# MW vs. LOGP
plot(NA, NA, xlim=c(0,600), ylim=c(-6.5,6.5), xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=(0:6)*100, lwd=1, lwd.ticks=1)
mtext(side=1, line=2.5, text='MW (Da)')
axis(side=2, at=-7:7, lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='log(P)')
rect(xleft=0, xright=300, ybottom=-6.5, ytop=3, border=NA, col=alpha(compliant_color, .5))
points(frags$mw, frags$logp, pch=20, cex=0.08, col='black')
compliant = mean(frags$mw <= 300 & frags$logp <= 3, na.rm=T)
text(x=150, y=-5, labels=paste(percent(compliant,digits=1)), col='black', font=2)
text(x=450, y=-5, labels=paste(percent(1-compliant,digits=1)), col='black', font=2)
mtext('B', side=3, cex=2, adj = -0.3, line = 0.5)

par(mar=c(3,4,3,1))

hba_hist = sqldf("select hba, count(*) n from frags group by 1 order by 1;")
hba_hist$dens = hba_hist$n / sum(hba_hist$n)
plot(NA, NA, xlim=c(0,10), ylim=c(0,.4), xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=c(0,10), labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, at=hba_hist$hba+0.5, labels=hba_hist$hba, lwd=0, lwd.ticks=0, cex.axis=0.8, line=-1)
mtext(side=1, line=1.5, text='HBA')
axis(side=2, at=(0:7)*5/100, labels=percent((0:7)*5/100,digits=0), lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=3, text='proportion')
compliant = hba_hist$hba <= 3
#rect(xleft=hba_hist$hba, xright=hba_hist$hba+1, ybottom=rep(0, nrow(hba_hist)), ytop=hba_hist$dens, border=NA, col='black')
rect(xleft=hba_hist$hba[!compliant], xright=hba_hist$hba[!compliant]+1, ybottom=rep(0, nrow(hba_hist[!compliant,])), ytop=hba_hist$dens[!compliant], border=NA, col='black')
rect(xleft=hba_hist$hba[compliant], xright=hba_hist$hba[compliant]+1, ybottom=rep(0, nrow(hba_hist[compliant,])), ytop=hba_hist$dens[compliant], border=NA, col=compliant_color)
compliant = mean(frags$hba <= 3, na.rm=T)
text(x=1.5, y=.325, labels=paste(percent(compliant,digits=1)), col='black', font=2)
text(x=6.0, y=.325, labels=paste(percent(1-compliant,digits=1)), col='black', font=2)
mtext('C', side=3, cex=2, adj = -0.3, line = 0.5)

hbd_hist = sqldf("select hbd, count(*) n from frags group by 1 order by 1;")
hbd_hist$dens = hbd_hist$n / sum(hbd_hist$n)
plot(NA, NA, xlim=c(0,10), ylim=c(0,.80), xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=c(0,10), labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, at=hba_hist$hba+0.5, labels=hba_hist$hba, lwd=0, lwd.ticks=0, cex.axis=0.8, line=-1) # still use HBA axis
mtext(side=1, line=1.5, text='HBD')
axis(side=2, at=(0:10)*10/100, labels=percent((0:10)*10/100,digits=0), lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=3, text='proportion')
compliant = hbd_hist$hbd <= 3
#rect(xleft=hbd_hist$hbd, xright=hbd_hist$hbd+1, ybottom=rep(0, nrow(hbd_hist)), ytop=hbd_hist$dens, border=NA, col='black')
rect(xleft=hbd_hist$hbd[!compliant], xright=hbd_hist$hbd[!compliant]+1, ybottom=rep(0, nrow(hbd_hist[!compliant,])), ytop=hbd_hist$dens[!compliant], border=NA, col='black')
rect(xleft=hbd_hist$hbd[compliant], xright=hbd_hist$hbd[compliant]+1, ybottom=rep(0, nrow(hbd_hist[compliant,])), ytop=hbd_hist$dens[compliant], border=NA, col=compliant_color)
compliant = mean(frags$hbd <= 3, na.rm=T)
text(x=1.5, y=.7, labels=paste(percent(compliant,digits=1)), col='black', font=2)
text(x=6.0, y=.7, labels=paste(percent(1-compliant,digits=1)), col='black', font=2)
mtext('D', side=3, cex=2, adj = -0.3, line = 0.5)

dev.off()
