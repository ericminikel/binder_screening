options(stringsAsFactors = F)
if(interactive()) {
  setwd('~/d/sci/src/binder_screening/')
}
library(sqldf)
library(rcdk)

normal_color = '#000123'
error_color = '#C9C9C9'
dmso_color = '#0276FD'

format_delta = function(x) {
  sign = ifelse(x == 0, '', ifelse(x > 0, '+', '-'))
  return( paste0(sign, abs(x)) )
}

percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='f'),"%",sep="") ) )
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot

tla_to_ola = function(x) {
  mapping = data.frame(tla=c("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","TER","THR","TRP","TYR","VAL"),
                       ola=c("A",  "R",  "N",  "D",  "C",  "Q",  "E",  "G",  "H",  "I",  "L",  "K",  "M",  "F",  "P",  "S",  "X",  "T",  "W",  "Y",  "V"))
  for (row in 1:dim(mapping)[1]) {
    x = gsub(mapping$tla[row], mapping$ola[row], toupper(x))
  }
  return (x)
}

#### TABLE 1 and FIGURE 1: FRAGMENT SCREENING

# parameters for table of libraries screened
liborder = data.frame(library=c('f19','chiral','cdot1','chembridge','cdot2'),
                      lorder=c(1,4,2,5,3),
                      disp=c('Broad 19F library','Schreiber chiral collection','Broad 1st gen STD library','Chembridge high-solubility set','Broad 2nd gen STD library'),
                      color=c("#66A61E", "#E7298A", "#7570B3", "#D95F02", "#1B9E77"))

frags = read.table('data/fragments/fragments.tsv', sep='\t', header=T, quote='', comment.char='')
frags$color = liborder$color[match(frags$library, liborder$library)]

# descriptors = get.desc.names(type="all")
properties = data.frame(
  descriptor = c("org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor",
                 "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor",             
                 "org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor",
                 "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor",
                 "org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor",
                 "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor"
  ),
  name = c("hba", "hbd", "logp", "mw", "rotb", "tpsa")
)

recompute_fragment_properties = F
if (recompute_fragment_properties) {
  for (i in 1:nrow(frags)) {
    write(paste0('\rComputing properties for compound ',i),stderr())
    flush.console()
    compound = parse.smiles(frags$smiles[i])[[1]]
    if (is.null(compound)) {
      frags$parsed[i] = FALSE
      next
    } else {
      frags$parsed[i] = TRUE
    }
    for (j in 1:nrow(properties)) {
      temp = eval.desc(compound, properties$descriptor[j])
      if (properties$descriptor[j] == "org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor") {
        temp = temp['ALogP']
      }
      frags[i, properties$name[j]] = round(temp,digits=1)
    }
  }
  write.table(frags, 'data/fragments/fragments_with_properties.tsv', sep='\t', col.names=T, row.names=F, quote=F)
} else {
  frags = read.table('data/fragments/fragments_with_properties.tsv', sep='\t', header=T, quote='', comment.char='')
}

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
               select   l.disp, count(*) n, count(distinct pool) hit_pools, sum(std_attempted) stds, sum(trosy_attempted) trosys, sum(trosy_validated) validated
               from     frags f, liborder l
               where    f.library = l.library
               group by 1
               order by l.lorder
               ;")

table1$n[table1$disp=='Schreiber chiral collection'] = table1$n[table1$disp=='Schreiber chiral collection'] + 40 # 40 undisclosable structures omitted from fragments.tsv
table1$stds[table1$disp=='Schreiber chiral collection'] = table1$stds[table1$disp=='Schreiber chiral collection'] + 1 # the 40 include 1 that was in the hit pool and went on to STD NMR
sum(table1$n)

write.table(table1, 'display_items/table1.tsv', sep='\t', row.names=F, col.names=T, quote=F)

table1$color = rev(c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e'))

table1$wrap = gsub('(.{1,20})(\\s|$)', '\\1\n', table1$disp) # https://stackoverflow.com/a/2352006/3806692
table1$wrap = trimws(table1$wrap) # remove terminal line breaks

resx=300
#### FIGURE 1 - properties of library
pdf('display_items/figure-1.pdf',width=6.5,height=2.75)

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













#### FIGURE 3: DSF SCREENING

dsf_primary = read.table('data/nibr_dsf/dsf_primary.tsv',sep='\t',header=T,comment.char='',quote='',as.is=T)

plate_stats1 = sqldf("select xp, file_number, median(tm_inf_) median_apo_tm from dsf_primary where xp is not null and dmso group by 1, 2 order by 1;")
plate_stats2 = sqldf("
select   p.xp, p.median_apo_tm, median(abs(d.tm_inf_ - p.median_apo_tm)) mad
from     dsf_primary d, plate_stats1 p
where    d.xp = p.xp
and      d.dmso
group by 1, 2
order by 1
;")

overall_median_tm = median(dsf_primary$tm_inf_[dsf_primary$dmso & !is.na(dsf_primary$xp)])
overall_mad = median(abs(dsf_primary$tm_inf_[dsf_primary$dmso & !is.na(dsf_primary$xp)] - overall_median_tm))
plate_stats2$delta_tm_median = plate_stats2$median_apo_tm - overall_median_tm

dsf_primary$delta_tm = dsf_primary$tm_inf_ - overall_median_tm
dsf_primary$delta_tm = pmax(dsf_primary$delta_tm,-1,pmin(dsf_primary$delta_tm,1))

dsf_primary$color = normal_color
dsf_primary$color[dsf_primary$error_code != 0] = error_color
dsf_primary$color[dsf_primary$dmso] = dmso_color

dsf_primary$delta_tm_disp = pmax(pmin(dsf_primary$delta_tm, 1),-1)

xlims = c(0,max(dsf_primary$xp, na.rm=T)+1)
xats = c(1, 20, 40, 60, 86)

dsf_retest = read.table('data/nibr_dsf/dsf_retest.tsv',sep='\t',header=T,quote='',comment.char='',as.is=T)
dsf_hsqc = read.table('data/nibr_dsf/dsf_hsqc.tsv',sep='\t',header=T,quote='',comment.char='',as.is=T)

dsf_retest$delta_tm = dsf_retest$tm_inf_ - dsf_retest$apo_tm_inf__median # vs. in-plate dmso controls, per Doug's email 2020-04-28

dsf_retest_avg = sqldf("
                       select   smiles, avg(delta_tm) mean_dtm, count(*), avg(delta_tm / apo_tm_inf__mad) mean_dtm_mads
                       from     dsf_retest
                       where    error_code = 0
                       and      not blank
                       and      smiles != ''
                       group by 1
                       order by 1
                       ;") # 1187
length(unique(dsf_retest_avg$smiles))
colnames(dsf_retest_avg)[2] = 'delta_tm'

dsf_primary$out3mad = (dsf_primary$tm_inf_ > overall_median_tm + 3*overall_mad | dsf_primary$tm_inf_ < overall_median_tm - 3*overall_mad)
dsf_primary_avg = sqldf("
                        select   smiles, avg(tm_inf_) mean_tm, count(*) n
                        from     dsf_primary
                        where    error_code = 0
                        and      picks != 'Untagged'
                        and      out3mad
                        and      smiles != ''
                        and      retested
                        group by 1
                        order by 1
                        ;") 
sum(dsf_primary_avg$mean_tm > overall_median_tm + 3*overall_mad)
sum(dsf_primary_avg$mean_tm < overall_median_tm - 3*overall_mad)

dsf_primary_avg$delta_tm = dsf_primary_avg$mean_tm - overall_median_tm

retest_results = sqldf("
                       select   p.delta_tm primary_delta, r.delta_tm retest_delta, p.smiles, r.mean_dtm_mads
                       from     dsf_primary_avg p, dsf_retest_avg r
                       where    p.smiles = r.smiles
                       ;") 
nrow(retest_results)
sum(retest_results$primary_delta < 0)
sum(retest_results$primary_delta > 0)

dsf_retest_neg_hits_smiles = read.table('data/nibr_dsf/dsf_retest_negative_hits.tsv',sep='\t',header=T,quote='',comment.char='',as.is=T)
dsf_retest_pos_hits_smiles = read.table('data/nibr_dsf/dsf_retest_positive_hits.tsv',sep='\t',header=T,quote='',comment.char='',as.is=T)

retest_results$retest_hit = 'no'
retest_results$retest_hit[retest_results$smiles %in% dsf_retest_neg_hits_smiles$smiles] = 'negative'
retest_results$retest_hit[retest_results$smiles %in% dsf_retest_pos_hits_smiles$smiles] = 'positive'
table(retest_results$retest_hit)
sum(retest_results$primary_delta < 0 & retest_results$mean_dtm_mads < -3)
sum(retest_results$primary_delta > 0 & retest_results$mean_dtm_mads > 3)
retest_results$primary_delta = pmax(pmin(retest_results$primary_delta, 3),-3)
retest_results$retest_delta = pmax(pmin(retest_results$retest_delta, 3),-3)

dsf_primary$retested = dsf_primary$smiles %in% retest_results$smiles & (dsf_primary$delta_tm > 3*overall_mad | dsf_primary$delta_tm < -3*overall_mad)
retest_results$hsqced = retest_results$smiles %in% dsf_hsqc$smiles
sum(retest_results$smiles %in% dsf_hsqc$smiles)

sum(retest_results$retest_delta[retest_results$hsqced] < 0)
sum(retest_results$retest_delta[retest_results$hsqced] > 0)


# identify the highlighted hit in the retest data
smiles_3c = dsf_primary$smiles[dsf_primary$file_number==54 & dsf_primary$wellnumber==139]

dsf_primary$hit3c = dsf_primary$smiles == smiles_3c
retest_results$hit3c = retest_results$smiles == smiles_3c

hit_color = '#FF2020'
resx=300


#png('display_items/figure-3a.png',width=6.5*resx,height=3*resx,res=resx)
cairo_pdf('display_items/figure-3a.pdf',width=6.5,height=2)
par(mar=c(3,4,3,4))
plot(NA, NA, xlim=xlims, ylim=c(-.5,.5), xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=xats)
axis(side=2, at=-5:5/10, labels=NA, las=2, tck=-0.025)
axis(side=2, at=-1:1/2, labels=format_delta(-1:1/2), las=2, tck=-0.05)
mtext(side=1, line=2.0, text='plate number')
mtext(side=2, line=2.5, text='ΔTm (°C)')
abline(h=0)
abline(h=c(-3,3)*overall_mad, lty=1, col=dmso_color)
mtext(side=4, at=c(-3,3)*overall_mad, line=0.25, text=paste0(c('-','+'),'3 MAD'), las=2, col=dmso_color, cex=0.8)
segments(x0=plate_stats2$xp-0.5,x1=plate_stats2$xp+0.5,y0=plate_stats2$delta_tm_median,col=dmso_color)
rect(xleft=plate_stats2$xp-0.5,xright=plate_stats2$xp+0.5,ybottom=plate_stats2$delta_tm_median-3*plate_stats2$mad,ytop=plate_stats2$delta_tm_median+3*plate_stats2$mad,col=alpha(dmso_color,ci_alpha), border=NA)
mtext('A', side=3, cex=2, adj = -0.1, line = 0.5)
dev.off()




#png('display_items/figure-3a.png',width=6.5*resx,height=3*resx,res=resx)
cairo_pdf('display_items/figure-3b.pdf',width=6.5,height=3)
par(mar=c(3,4,3,4))
plot(NA, NA, xlim=xlims, ylim=c(-1,1), xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=xats)
axis(side=2, at=-20:20/20, labels=NA, las=2, tck=-0.025)
axis(side=2, at=-2:2/2, labels=format_delta(-2:2/2), las=2, lwd=0, line=-0.5)
axis(side=2, at=-2:2/2, labels=NA, las=2, tck=-0.05)
mtext(side=1, line=2.0, text='plate number')
mtext(side=2, line=2.5, text='ΔTm (°C)')
abline(h=0)
points(dsf_primary$x[!dsf_primary$dmso], dsf_primary$delta_tm_disp[!dsf_primary$dmso], pch='.', cex=0.1, col=dsf_primary$color[!dsf_primary$dmso])
points(dsf_primary$x[dsf_primary$retested], dsf_primary$delta_tm_disp[dsf_primary$retested], pch='.', cex=0.5, col=hit_color)
par(xpd=T)
points(dsf_primary$x[dsf_primary$hit3c], dsf_primary$delta_tm_disp[dsf_primary$hit3c], pch=0, cex=1.5, col=hit_color)
par(xpd=F)
abline(h=c(-3,3)*overall_mad, lty=1, col=dmso_color)
mtext(side=4, at=c(-3,3)*overall_mad, line=0.25, text=paste0(c('-','+'),'3 MAD'), las=2, col=dmso_color, cex=0.8)
mtext('B', side=3, cex=2, adj = -0.1, line = 0.5)
dev.off()



#png('display_items/figure-3b.png',width=3.5*resx,height=3*resx,res=resx)
cairo_pdf('display_items/figure-3c.pdf',width=3.5,height=3)
par(mar=c(4,4,3,4))
plot(NA, NA, xlim=c(-3.1,3.1), ylim=c(-3.1,3.1), axes=F, ann=F, xaxt='n', yaxt='n', xaxs='i', yaxs='i')
axis(side=1, at=-3:3, labels=format_delta(-3:3), lwd=0, lwd.ticks=1, cex.axis=.75)
axis(side=2, at=-3:3, labels=format_delta(-3:3), lwd=0, lwd.ticks=1, las=2, cex.axis=.75)
#axis(side=1, at=c(-1,-.5,0,.5,1), lwd=0, lwd.ticks=1)
#axis(side=2, at=c(-1,-.5,0,.5,1), lwd=0, lwd.ticks=1, las=2)
mtext(side=1, line=2.5, text='primary screen ΔTm (°C)')
mtext(side=2, line=2.5, text='retest ΔTm (°C)')
abline(h=0)
abline(v=0)
#abline(h=c(-3,3)*overall_mad, col=dmso_color)
abline(v=c(-3,3)*overall_mad, col=dmso_color)
#mtext(side=4, line=0.25, at=c(-3,3)*overall_mad, col=dmso_color, las=2, text=c('-3 MAD','+3 MAD'), cex=0.75)
mtext(side=3, line=0.25, at=c(-3,3)*overall_mad, col=dmso_color, las=2, text=c('-3 MAD','+3 MAD'), cex=0.75)
abline(a=0, b=1, col='black')
points(retest_results$primary_delta, retest_results$retest_delta, col=normal_color, pch=20, cex=0.5)
points(retest_results$primary_delta[retest_results$hsqced], retest_results$retest_delta[retest_results$hsqced], col=hit_color, pch=20, cex=0.5)
par(xpd=T)
points(retest_results$primary_delta[retest_results$hit3c], retest_results$retest_delta[retest_results$hit3c], col=hit_color, pch=0, cex=2.0)
par(xpd=F)
mtext('C', side=3, cex=2, adj = -0.1, line = 0.5)
dev.off()


### DSF raw melt curves
# iconv -c -f utf-16 -t utf-8 spotfire_file > desired_file # https://stackoverflow.com/a/2398403/3806692
raw_primary = read.table('data/nibr_dsf/dsf_raw_primary.tsv',sep='\t',header=T,quote='',comment.char='',as.is=T)
raw_primary = raw_primary[raw_primary$file_number==54,]
hit_primary = raw_primary[raw_primary$well_number==139,]
dmso_primary = raw_primary[raw_primary$well_number %% 24 %in% c(0,23),] # this mod division gets cols 23-24 which are in-plate dmso controls

raw_retest = read.table('data/nibr_dsf/dsf_raw_retest.tsv',sep='\t',header=T,quote='',comment.char='',as.is=T)
dsf_retest[dsf_retest$smiles==smiles_3c,c('file_number','wellnumber')]
dsf_retest$file_well = paste0(dsf_retest$file_number,'-',dsf_retest$wellnumber)
hit_file_wells = dsf_retest[dsf_retest$smiles==smiles_3c,c('file_well')]
raw_retest$file_well = paste0(raw_retest$file_number,'-',raw_retest$well_number)

hit_retest = raw_retest[raw_retest$file_well %in% hit_file_wells,]

# also grab raw melt curves for in-plate dmso controls
# iconv -c -f utf-16 -t utf-8 data_received/HT.DSF_DMSO_curves.csv > data/dsf_retest_dmso_raw.csv
raw_retest_dmso = read.table('data/nibr_dsf/dsf_raw_retest_dmso.tsv',sep='\t',header=T,quote='',comment.char='',as.is=T)
colnames(raw_retest_dmso) = gsub('[^a-z0-9_]','_',tolower(colnames(raw_retest_dmso)))
raw_retest_dmso = raw_retest_dmso[raw_retest_dmso$file_number %in% hit_retest$file_number,]
mean(raw_retest_dmso$well_number %% 24 %in% c(0,23)) # check that all are final 2 columns - 1 is good
raw_retest_dmso$file_well = paste0(raw_retest_dmso$file_number,'-',raw_retest_dmso$well_number)


resx=300
#png('display_items/figure-3c.png',width=3.25*resx,height=3.25*resx,res=resx)
cairo_pdf('display_items/figure-3d.pdf',width=3.25,height=3.25)
par(mar=c(3,3,3,1))
plot(NA, NA, xlim=c(30,90), ylim=c(0,6.5), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=0:10*10)
axis(side=2, at=0:10, las=2)
mtext(side=1, line=2, text='temperature (°C)')
mtext(side=2, line=2, text='intensity')
for (w in unique(dmso_primary$well_number)) {
  subs = subset(dmso_primary, well_number == w)
  points(x=subs$temperature, y=subs$intensity, type='l', lwd=0.125, col=alpha(dmso_color,.1))
}
for (fw in unique(raw_retest_dmso$file_well)) {
  subs = subset(raw_retest_dmso, file_well == fw)
  points(x=subs$temperature, y=subs$intensity, type='l', lwd=0.125, col=alpha(dmso_color,.1))
}
dmso_all = rbind(dmso_primary[,c('temperature','intensity')], raw_retest_dmso[,c('temperature','intensity')])
dmso_meancurve = sqldf("select temperature, avg(intensity) mean_intens, count(*) n from dmso_all group by 1 order by 1;")
points(x=dmso_meancurve$temperature, y=dmso_meancurve$mean_intens, type='l', lwd=2, col=dmso_color)
points(x=hit_primary$temperature, y=hit_primary$intensity, type='l', lwd=0.25, col=alpha(hit_color,.4))
for (w in unique(hit_retest$file_well)) {
  subs = subset(hit_retest, file_well == w)
  points(x=subs$temperature, y=subs$intensity, type='l', lwd=0.25, col=alpha(hit_color,.4))
}
hit_all = rbind(hit_primary[,1:5], hit_retest[,1:5])
hit_meancurve = sqldf("select temperature, avg(intensity) mean_intens, count(*) n from hit_all group by 1 order by 1;")
points(x=hit_meancurve$temperature, y=hit_meancurve$mean_intens, type='l', lwd=2, col=hit_color)

dmso_ymed_initial = (max(dmso_meancurve$mean_intens) + min(dmso_meancurve$mean_intens))/2
dmso_which = which(dmso_meancurve$mean_intens >= dmso_ymed_initial)[1]
dmso_tmed = approx(y=dmso_meancurve$temperature[c(dmso_which-1, dmso_which)], x=dmso_meancurve$mean_intens[c(dmso_which-1, dmso_which)], xout=dmso_ymed_initial)$y
dmso_ymed = approx(x=dmso_meancurve$temperature[c(dmso_which-1, dmso_which)], y=dmso_meancurve$mean_intens[c(dmso_which-1, dmso_which)], xout=dmso_tmed)$y
hit_ymed_initial = (max(hit_meancurve$mean_intens) + min(hit_meancurve$mean_intens))/2
hit_which = which(hit_meancurve$mean_intens >= hit_ymed_initial)[1]
hit_tmed = approx(y=hit_meancurve$temperature[c(hit_which-1, hit_which)], x=hit_meancurve$mean_intens[c(hit_which-1, hit_which)], xout=dmso_ymed_initial)$y
mean_ymed = (dmso_ymed_initial + hit_ymed_initial)/2

arrows(x0=dmso_tmed, x1=hit_tmed, y0=mean_ymed, code=3, angle=90, length=0.025, col='black')
delta_tm = round(hit_tmed - dmso_tmed,1)
text(x=hit_tmed, y=dmso_ymed-0.5, pos=4, labels=paste0('ΔTm\n',delta_tm,'°C'), cex=0.8)

leg = data.frame(disp=c('DMSO replicate','DMSO mean','hit replicate','hit mean'),
                 col=rep(c(dmso_color,hit_color),each=2),
                 lwd=rep(c(.25,2),2))
legend('topleft',legend=leg$disp,col=leg$col,lwd=leg$lwd,text.col=leg$col,bty='n',cex=0.8)
mtext('D', side=3, cex=2, adj = -0.1, line = 0.5)
dev.off()

# check final stats
length(unique(dsf_primary$smiles))
nrow(retest_results)
sum(retest_results$hsqced)












resx=600

del = read.table('data/del/macrocycle_enrichment.tsv',sep='\t',header=T,quote='',comment.char='')

# this one has to be raster because so many points, vector becomes too slow to load
png('display_items/figure-4a.png',width=6.5*resx,height=3*resx,res=resx)
par(mar=c(4,4,3,1))

default_color = '#000000'
x_krd_color = '#BBBBBB'
x_cjs_color = '#00B050'
xcc_s_color = '#00B0F0'

leg = data.frame(leg=c('all compounds','*CJS series','CC*S series','*KRD series'),
                 col=c(default_color,x_cjs_color,xcc_s_color,x_krd_color))

del$prerank = rank(del$pre_enrich_count) # this ranks 0 counts as highest, max counts as lowest - but axis will be reversed below
del$enrich = del$raw_selection_freq/del$pre_enrich_freq
ymax = 250
ylims = c(0,ymax*1.05)
del$enrich[del$enrich > ymax] = ymax
del$color = default_color
#del$color[del$sequence %in% c('OCJS','NCJS')] = x_cjs_color
del$color[grepl('[A-Z]CJS',del$sequence)] = x_cjs_color
#del$color[del$sequence %in% c('CCHS','CCSS','CCTS','CCLS')] = xcc_s_color
del$color[grepl('CC[A-Z]S',del$sequence)] = xcc_s_color
del$color[grepl('[A-Z]KRD',del$sequence)] = x_krd_color
to_label = del$sequence %in% c('OCJS','NCJS','CCHS','CCSS','CCTS','CCLS','HKRD','FKRD')
xlims = c(0,max(del$prerank))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
xats = c(1,seq(5e4,2e5,by=5e4),256000)
axis(side=1, at=xats, labels=NA)
axis(side=1, at=xats, labels=formatC(rev(xats),format='d',big.mark=','), lwd=0, line=-0.5) # note reverse axis - lowest abundance at left, highest at right
mtext(side=1, line=1.5, text='pre-enrichment rank abundance')
axis(side=2, at=0:5*50, labels=c('0','50','100','150','200','250+'),las=2)
mtext(side=2, line=3.0, text='fold enrichment')
points(del$prerank, del$enrich, col=del$color, pch=20, cex=0.2)
points(del$prerank[del$color != default_color], del$enrich[del$color != default_color], col=del$color[del$color != default_color], pch=20, cex=1)
text(del$prerank[to_label], del$enrich[to_label], labels=del$sequence[to_label], pos=2, col=del$color[to_label], cex=0.9)
par(xpd=T)
legend(x=0,y=-65,legend=leg$leg,col=leg$col,text.col=leg$col,pch=19,horiz=T,bty='n',cex=0.9)
par(xpd=F)

mtext(side=3,adj=0,padj=-0.5,line=0.5,font=2,cex=2,text='A')
dev.off()






awise_dsf = read.table('data/awise/dsf.tsv',sep='\t',header=T)
awise_dsf$sem = awise_dsf$sd_delta_tm/sqrt(awise_dsf$n)
awise_dsf$l95 = awise_dsf$average_delta_tm - 1.96*awise_dsf$sem
awise_dsf$u95 = awise_dsf$average_delta_tm + 1.96*awise_dsf$sem

awise_dsf$x = awise_dsf$compound

awise_apo_sd = 0.263 # from Andrew's Excel file 2020-04-24
cutoffs = c(-3, 3)*awise_apo_sd

default_color = '#777777'
hit_color = '#FF2020'

xlims = range(awise_dsf$x) + c(-0.5, 0.5)
ylims = c(-1, 1)
awise_dsf$pch = 19 # default is a normal point
awise_dsf$pch[awise_dsf$average_delta_tm > max(ylims)] = 24
awise_dsf$pch[awise_dsf$average_delta_tm < min(ylims)] = 25
awise_dsf$dtm_disp = pmax(pmin(awise_dsf$average_delta_tm,1),-1)
awise_dsf$col = default_color

cairo_pdf('display_items/figure-5a.pdf',width=6.5,height=1.75)
par(mar=c(3,3,2.5,3))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=1:nrow(awise_dsf), labels=NA, tck=-0.025)
axis(side=1, at=seq(1,nrow(awise_dsf),by=20), labels=NA, tck=-0.05)
axis(side=1, at=seq(1,nrow(awise_dsf),by=20), lwd=0, line=-0.5)
axis(side=2, at=seq(-1,1,.5), labels=NA)
axis(side=2, at=seq(-1,1,.5), line=-0.5, las=2, lwd=0)
abline(h=0)
abline(h=cutoffs, col=dmso_color)
mtext(side=1, line=1.75, text='compound number')
mtext(side=2, line=2, text='ΔTm (°C)')
segments(x0=awise_dsf$x, y0=awise_dsf$l95, y1=awise_dsf$u95, lwd=2, col=awise_dsf$col)
par(xpd=T)
points(x=awise_dsf$x, y=awise_dsf$dtm_disp, pch=awise_dsf$pch, col=awise_dsf$col, bg=awise_dsf$col)
par(xpd=F)
mtext(side=4, line=0.25, at=cutoffs, col=dmso_color, text=paste0(c('-','+'),'3 SD'), las=2)
mtext('A', side=3, cex=2, adj = -0.1, line = 0.5)
dev.off()




#### FIGURE S1: TROSY characterization of HuPrP90-231

resx=600
png('display_items/figure-s1.png',width=6.5*resx,height=4*resx,res=resx)

# first run: src/bruker_2d_data.py -p data/nmr/190507_HuPrP90-231_TROSY_DMSO.txt > data/nmr/190507_HuPrP90-231_TROSY_DMSO.matrix
path = 'data/nmr/190507_HuPrP90-231_TROSY_DMSO.matrix'

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

calzolai = read_bmrb_hsqc('data/nmr/BMRB_5713.txt',add=118)

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
# here use expression() to get superscript ( https://stackoverflow.com/a/10628694/3806692 )
# because default PNG renderer recognizes UTF-8 superscript 1 but not superscript 5:
mtext(side=1, line=2.5, text=expression(' '^'1'*'H chemical shift (ppm)'))
mtext(side=2, line=3.5, text=expression(' '^'15'*'N chemical shift (ppm)'))
par(xpd=T)
text(x=calzolai$h1+.2, y=calzolai$n15-1, labels=calzolai$label, col=text_color, pos=3, font=2, cex=0.6)
par(xpd=F)

dev.off()






#### FIGURE S4B: AKTA ELUTION CURVES

pdf('display_items/figure-s4b.pdf',width=3.5,height=3.5)

relpath = 'data/elutions/'
elution_files = list.files(relpath)
flowrate = 6 # mL/min flow rate during AKTA elutions

timemax = 20
maumax = 1500
concmax = .5
par(mar=c(3,4,3,4))
plot(NA, NA, xlim=c(0,timemax), ylim=c(0,maumax), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=0:timemax, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(0,timemax,by=5), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=seq(0,timemax,by=5), labels=seq(0,timemax,by=5), lwd=0, line=-0.5)
mtext(side=1, line=1.75, text='time (minutes)', font=1)
axis(side=2, at=seq(0,maumax,by=100), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=2, at=seq(0,maumax,by=500), lwd=0, lwd.ticks=1, tck=-0.05, las=2)
mtext(side=2, line=3.0, text='A280 (mAU)', col='#0001CD', font=1)
# 
for (elution_file in elution_files) {
  path = paste0(relpath,elution_file)
  elution_raw = read.table(path,sep='\t',header=T,quote='"',comment.char='',fill=T)
  
  elution_uv = elution_raw[4:nrow(elution_raw),1:2]
  colnames(elution_uv) = tolower(elution_raw[2,1:2])
  elution_uv$min = as.numeric(elution_uv$min)
  elution_uv$mau = as.numeric(elution_uv$mau)
  
  elution_concb = elution_raw[4:nrow(elution_raw),5:6]
  colnames(elution_concb) = tolower(elution_raw[2,5:6])
  elution_concb$min = as.numeric(elution_concb$min)
  elution_concb$percentb = as.numeric(elution_concb[,'%'])/100
  
  elution_frac = elution_raw[3:9,9:10]
  colnames(elution_frac) = tolower(elution_raw[2,9:10])
  elution_frac$min = as.numeric(elution_frac$min)
  
  points(elution_uv$min, elution_uv$mau, col='#0001CD', type='l', lwd=1)

}

par(new=T)
plot(NA, NA, xlim=c(0,timemax), ylim=c(0,concmax), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=4, at=seq(0,concmax,by=.1), labels=percent(seq(0,concmax,by=.1),digits=0), lwd=1, lwd.ticks=1, las=2)
mtext(side=4, line=3.0, text='imidazole gradient (%B)', col='#FF9912', font=1)
# the conc b curve is same for all six, so no need to loop - just plot the final one:
points(elution_concb$min, elution_concb$percentb, col='#FF9912', type='l', lwd=2)
mtext('B', side=3, cex=2, adj = -0.1, line = 1.5)

dev.off()






