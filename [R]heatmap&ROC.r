require(tidyverse)
require(pheatmap)
require(matrixStats)

make.rownames=function (df, rowname.column = 1, keep.rowname.column = T) 
{
    df = as.data.frame(df)
    rownames(df) = df[[rowname.column]]
    if (!keep.rowname.column) {
        df[[rowname.column]] = NULL
    }
    df
}
extract.between=function (x, leading = "^", trailing = "$") 
{
    sub(paste0(pattern = ".*", leading, " *(.*?) *", trailing, 
        ".*"), replacement = "\\1", x)
}
mat = read.csv('/NAS/wg_zsx/notebooks/projects/wanghy/3sequence_IDH+TERT_features.csv',header = T)[-2] %>% make.rownames('id',keep.rowname.column = F)

feature_names = mat %>% colnames

feature_annot = data.frame(row.names = feature_names,
                           transformation=feature_names %>% extract.between('_','_')%>%extract.between('^','\\.'),
                           sequence = feature_names %>% extract.between('^','_'),
                           type = case_when(
                               grepl('shape',feature_names)~'shape',
                               grepl('firstorder',feature_names)~'firstorder',
                               TRUE~'texture'
                           ))

colData = read.csv('/NAS/wg_zsx/notebooks/projects/wanghy/IDH+TERT_heatmap.csv',header = T)[1:4] %>% make.rownames('id',keep.rowname.column = F)

colData = colData %>% within({
    idh = ifelse(idh==1,'IDHmut','IDHwt')
    tert = ifelse(tert==1,'pTERTmut','pTERTwt')
    idh.tert = ifelse(idh.tert==1,'double','others')
})

palette = list(idh=c(IDHmut='gold',IDHwt='white'),tert=c(pTERTmut='purple',pTERTwt='white'),
               idh.tert=c(double='royalblue','others'='white'),
               sequence=c(t1Gd='red',flair='skyblue',ADC='forestgreen'),
               transformation = c(original = 'skyblue',wavelet='orange',log='gold')
              )

mat_order=function (mat, sort = T,row.dist=col_cor_dist,col.dist=col_dist,hc.method='ward.D2') 
{
    mat[var_order(mat, sort,row.dist,hc.method), sample_order(mat, sort,col.dist,hc.method)]
}

col_cor_dist = function(mat){
    1 - cor(mat %>% t) %>% as.dist
}

col_dist = function(mat){(mat %>% t) %>% dist}

var_order=function (mat, sort = T,fun.dist=col_cor_dist,hc.method = 'ward.D2') 
{
    md = fun.dist(mat)
    hc = hclust(md, hc.method)
    if (sort) {
        hc = hc %>% dendsort
    }
    hc %>% hc_ordered
}

sample_order=function (mat, sort = T,fun.dist=cor_dist,hc.method = 'ward.D2') 
{
     md = fun.dist(mat)
    hc = hclust(md, hc.method)
    if (sort) {
        hc = hc %>% dendsort
    }
    hc %>% hc_ordered
}

rowZscore = function (x) 
{
    x = as.matrix(x)
    (x - rowMeans(x))/rowSds(x)
}

outlier.trimmer4heatmap=function (x, outlier.threshold = 1.5) 
{
    x = as.matrix(x)
    rowMeidan = rowMedians(x)
    offset = outlier.threshold * rowIQRs(x)
    top = matrix(rowMeidan + offset, nrow = nrow(x), ncol = ncol(x))
    bottom = matrix(rowMeidan - offset, nrow = nrow(x), ncol = ncol(x))
    x = ifelse(x > top, top, x)
    x = ifelse(x < bottom, bottom, x)
    x
}

htmap = mat %>% t %>% rowZscore %>% mat_order(hc.method = 'ward.D2') %>% outlier.trimmer4heatmap(1.5)%>%
  pheatmap(cluster_cols = F, cluster_rows = F, border_color = NA,width = 10,height=6,filename = '/NAS/wg_zsx/notebooks/projects/wanghy/htmap.pdf',
    silent = T, legend = F, #annotation_legend = F, 
    annotation_col = colData[c(3:1)], annotation_row = feature_annot,
    show_rownames = F, show_colnames = F, 
    annotation_colors = palette, 
   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100))

options(repr.plot.width = 10, repr.plot.height=6)
htmap

cormap = mat %>% cor(method = 'spearman') %>% mat_order(hc.method = 'ward.D2') %>% 
  pheatmap(cluster_cols = F, cluster_rows = F, border_color = NA,width = 5,height=5,filename = '/NAS/wg_zsx/notebooks/projects/wanghy/cormap.pdf',
    silent = T, legend = F, #annotation_legend = F, 
    annotation_col = feature_annot, annotation_row = feature_annot,
    show_rownames = F, show_colnames = F, 
    annotation_colors = palette, 
   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100))

options(repr.plot.width = 5, repr.plot.height=5)
cormap


require(pROC)

require(tidyverse)

require(MLmetrics)
require(ROCR)

df.split=function (df, f, sep = "-") 
{
    if (length(f) > 1) {
        split.f = df[, f, drop = F] %>% Reduce(f = function(x, 
            y) {
            paste0(x, sep, y)
        })
    }
    else {
        split.f = df[[f]]
    }
    df %>% split(f = split.f)
}

sequences = c(t1='red',flair='skyblue',ADC='forestgreen')

multimodal = read.csv('/NAS/wg_zsx/notebooks/projects/wanghy/multimodal_auc.csv')%>% df.split('s')

multimodal.stats = multimodal %>% names %>% flexible.mclapply(mc.cores = 200,function(m){
    with(multimodal[[m]],{
        data.frame(
            stringsAsFactors = F,
            resample = m,
            PRAUC=PRAUC(p, target),
            AUC=AUC(p, target),
            F1_Score=F1_Score(as.character(target), as.character(y), positive = '1'),
            Precision=Precision(as.character(target), as.character(y), positive = '1'),
            Recall=Recall(as.character(target), as.character(y), positive = '1'),
            Specificity = Specificity(as.character(target), as.character(y), positive = '1'),
            Sensitivity = Sensitivity(as.character(target), as.character(y), positive = '1'),
            Accuracy = Accuracy(as.character(target), as.character(y))
        )
    })
})

multimodal.stats = multimodal.stats[sapply(multimodal.stats,class)!='try-error']

multimodal.stats = multimodal.stats %>% bind_rows

multimodal.preds = multimodal[multimodal.stats$resample] %>% lapply(function(d){
    d %>% with(ROCR::prediction(p, target))
})

options(repr.plot.width=5,repr.plot.height=5)
plot(performance(multimodal.preds[['125']], "tpr", "fpr"))
text(x = 0.8,y=0.2,sprintf("AUC = %.3f",multimodal.stats[multimodal.stats$resample=='125','AUC']))

options(repr.plot.width=5,repr.plot.height=5)
plot(performance(multimodal.preds[['125']], "prec", "rec"))
text(x = 0.4,y=0.3,sprintf("PRAUC = %.3f",multimodal.stats[multimodal.stats$resample=='125','PRAUC']))

options(repr.plot.width=5,repr.plot.height=5)
pdf('/NAS/wg_zsx/notebooks/projects/wanghy/roc-pr.pdf')
plot(performance(multimodal.preds[['125']], "tpr", "fpr"))
text(x = 0.8,y=0.2,sprintf("AUC = %.3f",multimodal.stats[multimodal.stats$resample=='125','AUC']))
plot(performance(multimodal.preds[['125']], "prec", "rec"))
text(x = 0.4,y=0.3,sprintf("PRAUC = %.3f",multimodal.stats[multimodal.stats$resample=='125','PRAUC']))
dev.off()

require(ggpubr)

multimodal.imp = read.csv('/NAS/wg_zsx/notebooks/projects/wanghy/multimodal_imp.csv',stringsAsFactors = F)%>% within({
    feature = factor(feature,sapply(split(weight,feature),mean)%>%sort()%>%names)
})

options(repr.plot.width=15,repr.plot.height=8)
multimodal.imp %>% ggbarplot(x = 'feature',y='weight',add = "mean_se") + coord_flip()

pdf('/NAS/wg_zsx/notebooks/projects/wanghy/imp.pdf',width = 15,height = 8)
multimodal.imp %>% ggbarplot(x = 'feature',y='weight',add = "mean_se") + coord_flip()
dev.off()


