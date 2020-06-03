options(stringsAsFactors = F)
if (interactive()) {
  setwd('~/d/sci/src/fragments_prp')
}

# imgmode = 'png'
imgmode = 'pdf'
imgsave = get(imgmode) # get the function that saves as that type of image
if (imgmode == 'png') {
  resx = 600 # multiplier for image width, height, and resolution
} else if (imgmode == 'pdf') {
  resx = 1
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='f'),"%",sep="") ) )
}

tla_to_ola = function(x) {
  mapping = data.frame(tla=c("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","TER","THR","TRP","TYR","VAL"),
                       ola=c("A",  "R",  "N",  "D",  "C",  "Q",  "E",  "G",  "H",  "I",  "L",  "K",  "M",  "F",  "P",  "S",  "X",  "T",  "W",  "Y",  "V"))
  for (row in 1:dim(mapping)[1]) {
    x = gsub(mapping$tla[row], mapping$ola[row], toupper(x))
  }
  return (x)
}

imgsave(paste('figures/figure-2.',imgmode,sep=''),width=6.5*resx,height=6.5*resx)

layout_matrix = matrix(c(1,2,3,3),nrow=2,byrow=T)
layout(layout_matrix, heights=c(1,1.75))

elution_raw = read.table('data/190429_HuPrP90-231_prep67_15N_elute_001.tsv',sep='\t',header=T,quote='"',comment.char='')

flowrate = 6 # ml/min

elution_uv = elution_raw[4:nrow(elution_raw),1:2]
colnames(elution_uv) = tolower(elution_raw[2,1:2])
elution_uv$ml = as.numeric(elution_uv$ml)
elution_uv$mau = as.numeric(elution_uv$mau)
elution_uv$time = elution_uv$ml / flowrate

elution_concb = elution_raw[4:nrow(elution_raw),5:6]
colnames(elution_concb) = tolower(elution_raw[2,5:6])
elution_concb$ml = as.numeric(elution_concb$ml)
elution_concb$percentb = as.numeric(elution_concb[,'%'])/100
elution_concb$time = elution_concb$ml / flowrate

elution_frac = elution_raw[3:9,9:10]
colnames(elution_frac) = tolower(elution_raw[2,9:10])
elution_frac$ml = as.numeric(elution_frac$ml)
elution_frac$time = elution_frac$ml / flowrate

timemax = 25
maumax = 1000
concmax = .5
par(mar=c(4,5,4,5))
plot(NA, NA, xlim=c(0,timemax), ylim=c(0,maumax), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=0:timemax, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(0,timemax,by=5), labels=seq(0,timemax,by=5), lwd=0, lwd.ticks=1, tck=-0.05)
mtext(side=1, line=2.5, text='time (minutes)', font=1)
axis(side=2, at=seq(0,maumax,by=100), lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=3.0, text='A280 (mAU)', col='#0001CD', font=1)
points(elution_uv$time, elution_uv$mau, col='#0001CD', type='l', lwd=2)
par(new=T)
plot(NA, NA, xlim=c(0,timemax), ylim=c(0,concmax), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=4, at=seq(0,concmax,by=.1), labels=percent(seq(0,concmax,by=.1),digits=0), lwd=1, lwd.ticks=1, las=2)
mtext(side=4, line=3.0, text='imidazole gradient (%B)', col='#FF9912', font=1)
points(elution_concb$time, elution_uv$percentb, col='#FF9912', type='l', lwd=1)
keptfrac = which(elution_frac$fraction=='1.A.3')
abline(v=elution_frac$time[keptfrac:(keptfrac+1)],col='red',lwd=0.5)
mtext(side=3, at=mean(elution_frac$time[keptfrac:(keptfrac+1)]), line=0.5, col='red', text='fraction\ncollected')
mtext('A', side=3, cex=2, adj = -0.1, line = 1.5)

# blank B for gel
par(mar=c(4,5,4,5))
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), axes=F, ann=F, xaxs='i', yaxs='i')
mtext('B', side=3, cex=2, adj = -0.1, line = 1.5)


par(mar=c(4,5,4,5))
# first run src/bruker_2d_data.py -p data/190507_HuPrP90-231_TROSY_DMSO.txt > data/190507_HuPrP90-231_TROSY_DMSO.matrix
path = 'data/190507_HuPrP90-231_TROSY_DMSO.matrix'

spec = as.matrix(read.table(path,row.names=1,header=TRUE,quote='',comment.char='',skip=8,sep='\t'))
spec = pmax(spec, 0.0)
spec = pmin(spec, 20000)
spec = t(spec)
# check what the matrix looks like
# spec[1:3,1:3]

n15vals = -as.numeric(gsub('y','',colnames(spec)))
h1vals = -as.numeric(gsub('x','',rownames(spec)))

read_bmrb_hsqc = function(path, add=0) {
  bmrb = read.table(path, header=FALSE)
  hsqc = merge(bmrb[bmrb$V8=='N',c('V6','V7','V11')], bmrb[bmrb$V8=='H',c('V6','V11')], by='V6')
  colnames(hsqc) = c('num','resi','n15','h1')
  hsqc$n15 = -hsqc$n15
  hsqc$h1 = -hsqc$h1
  hsqc$aa = tla_to_ola(hsqc$resi)
  hsqc$aanum = add + as.integer(hsqc$num)
  hsqc$label = paste(hsqc$aa, hsqc$aanum, sep='')
  return (hsqc)
}

calzolai = read_bmrb_hsqc('data/BMRB_5713.txt',add=118)

line_color = '#777777'
axis_color = '#000000'
spec_color = '#0001CD'
text_color = "#FF2019"

par(mar=c(4,5,3,1))
plot(NA, NA, axes=F, ann=F, xlim=c(-10.3, -6.2), ylim=c(-128, -107), xaxs='i', yaxs='i')
xats = -11:-4
yats = seq(-105,-135,by=-5)
abline(h=yats, col=line_color)
abline(v=xats, col=line_color)
axis(side=1, at=xats, lwd=0, lwd.ticks=0 ,col.axis=axis_color)
axis(side=2, at=yats, lwd=0, lwd.ticks=0, las=2, col.axis=axis_color)
levels = (6:7)*1000
contour(x=h1vals, y=n15vals, spec, levels=levels, axes=FALSE, ann = FALSE, add=TRUE, col=alpha(spec_color,.35))
levels = (8:20)*1000
contour(x=h1vals, y=n15vals, spec, levels=levels, axes=FALSE, ann = FALSE, add=TRUE, col=spec_color)
mtext(side=1, line=2.5, text='1H chemical shift (ppm)')
mtext(side=2, line=3.5, text='15N chemical shift (ppm)')  
par(xpd=T)
text(x=calzolai$h1+.2, y=calzolai$n15-1, labels=calzolai$label, col=text_color, pos=3, font=2, cex=0.6)
par(xpd=F)
mtext('C', side=3, cex=2, adj = 0.0, line = 0.5)

dev.off()