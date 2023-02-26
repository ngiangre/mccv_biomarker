# Setup -------------------------------------------------------------------

pacman::p_load(tidyverse,data.table)

mytheme <-
    ggthemes::theme_clean(base_size = 16) +
    theme(
        legend.position = "none"
    )

img_dir="docs/imgs/"

# Load datasets `datasets` -----------------------------------------------------------

dists = c("normal","t","beta")
datasets = lapply(dists,function(dist){
    tmp = fread(paste0("data/mccv_50subjects_500parameters_",dist,"_data.csv")) %>%
        mutate(class=factor(class,levels=c(0,1)))
    colnames(tmp)[1] = "subject"
    tmp %>%
        .[order(seed,subject)]
})
names(datasets) <- dists


# Compute ttests and FC statistics on datasets ---------------------------------

datasets_summary <- lapply(dists,function(dist_){
    seeds=datasets[[dist_]][,unique(seed)]
    map(seeds,function(seed_){
        a=datasets[[dist_]][seed==seed_ & class==0,biomarker]
        b=datasets[[dist_]][seed==seed_ & class==1,biomarker]
        t_test = t.test(b,a)
        t_stat = t_test$statistic
        t_pval = t_test$p.value
        data.table(seed=seed_,class0_mean=mean(a),class1_mean=mean(b),
                   t_statistic=t_stat,t_pvalue=t_pval,
                   fc=mean(b)/mean(a),
                   fc_median=median(b)/median(a),
                   fc_maxmin=min(b)/max(a),
                   log2fc=log2(mean(b)) - log2(mean(a)) )
    }) %>% bind_rows()
})
names(datasets_summary) <- dists


# Volcano plot of distribution tests --------------------------------------

g <-
    datasets_summary %>%
    bind_rows(.id="distribution") %>%
    arrange(t_statistic) %>%
    ggplot(aes(t_statistic,-log10(t_pvalue))) +
    geom_point(pch=21,size=3) +
    facet_wrap(~distribution,ncol=1) +
    colorspace::scale_fill_continuous_diverging(mid=0.75) +
    xlab("T statistic") +
    ylab("-log10(pvalue)") +
    ggthemes::theme_clean()
ggsave(paste0(img_dir,"ttest_statistic_vs_significance_by_distributions.pdf"),
       width=6,height=9)

# Load MCCV data `mccv_data` ----------------------------------------------------------

dists = c("normal","t","beta")
mccv_data = lapply(dists,function(dist){
    dir_ = paste0("data/mccv_output/",dist,"/real/")
    files = list.files(dir_)
    filenames = map(files,str_split,pattern="\\.") %>% map_chr(c(1,1))
    filetypes = map(filenames,str_split,pattern="_",n=3) %>% map_chr(c(1,3))
    seeds = map(filenames,str_split,pattern="data") %>% map_chr(c(1,1))
    lst = list()
    for(ft in unique(filetypes)){
        lst[[ft]] = lapply(files[filetypes==ft],
                           function(x){fread(paste0(dir_,x),nThread = 4)}
                           ) %>% bind_rows()
    }
    lst
})
names(mccv_data) <- dists



# Load MCCV permuted data `mccv_permuted_data` ----------------------------------------------------------

dists = c("normal","t","beta")
mccv_permuted_data = lapply(dists,function(dist){
    dir_ = paste0("data/mccv_output/",dist,"/permuted/")
    files = list.files(dir_)
    filenames = map(files,str_split,pattern="\\.") %>% map_chr(c(1,1))
    filetypes = map(filenames,str_split,pattern="_",n=4) %>% map_chr(c(1,4))
    seeds = map(filenames,str_split,pattern="data") %>% map_chr(c(1,1))
    lst = list()
    for(ft in unique(filetypes)){
        lst[[ft]] = lapply(files[filetypes==ft],
                           function(x){fread(paste0(dir_,x),nThread = 4)}
        ) %>% bind_rows()
    }
    lst
})
names(mccv_permuted_data) <- dists

# `mccv_biomarker_boxplot` and  `mccv_prediction_boxplot` -------------------------------------------------

mccv_biomarker_boxplot <- function(i,dist="normal",theme=mytheme){
    datasets[[dist]][seed==i] %>%
        ggplot(aes(class,biomarker,fill=class,group=class)) +
        geom_boxplot(outlier.shape = NA,fill="gray95") +
        geom_point(pch=21,size=3,
                   position = ggbeeswarm::position_quasirandom()) +
        scale_fill_brewer(palette = "Set1") +
        xlab("Class") +
        ylab("Biomarker expression") +
        theme
}

mccv_prediction_boxplot <- function(i,dist="normal",modelname=c("Logistic Regression"),theme=mytheme){
    mccv_data[[dist]][["patient_predictions"]][seed==i] %>%
        filter(model==modelname) %>%
        rename(subject=`V1`) %>%
        .[,.(lwr=quantile(y_proba,0.025,na.rm=TRUE),
             mean=mean(y_proba,na.rm=TRUE),
             upr=quantile(y_proba,0.025,na.rm=TRUE)),
          .(model,seed,subject,y_true)] %>%
        mutate(y_true = factor(y_true,levels=c(0,1))) %>%
        ggplot(aes(y_true,mean,fill=y_true,group=y_true)) +
        geom_boxplot(outlier.shape = NA,fill="gray95") +
        geom_point(pch=21,size=3,
                   position = ggbeeswarm::position_quasirandom()) +
        scale_fill_brewer(palette = "Set1") +
        xlab("Class") +
        ylab("Subject prediction") +
        ylim(0,1) +
        theme
}


# `mccv_performance_summary`, `mccv_prediction_summary`, `mccv_importance_summary` ------------------------------------------------------

mccv_performance_summary <- lapply(dists,function(dist){
    mccv_data[[dist]]$performance[,.(lwr = quantile(value,0.025,na.rm=TRUE),
                                            mean = mean(value,na.rm=TRUE),
                                            upr = quantile(value,0.975,na.rm=TRUE)),
                                         .(model,metric,seed)] %>%
        .[order(mean)] %>%
        mutate(significant = lwr>0.5 | upr<0.5)
})
names(mccv_performance_summary) <- dists

mccv_prediction_summary <- lapply(dists,function(dist_){
    mccv_data[[dist_]]$patient_prediction[,.(
        class0_mean = mean(y_proba[y_pred==0]),
        class1_mean = mean(y_proba[y_pred==1]),
        fc = mean(y_proba[y_pred==1])/mean(y_proba[y_pred==0]),
        t_statistic = t.test(y_proba[y_pred==1],
                          y_proba[y_pred==0])[1],
        t_pvalue = t.test(y_proba[y_pred==1],
                        y_proba[y_pred==0])[3]
    ),
        .(model,seed)] %>%
        .[order(fc,decreasing = TRUE)] %>%
        mutate(significant = fc>2)
})
names(mccv_prediction_summary) <- dists

mccv_importance_summary <- lapply(dists,function(dist){
    mccv_data[[dist]]$feature_importance[,.(lwr = quantile(importance,0.025,na.rm=TRUE),
                                           mean = mean(importance,na.rm=TRUE),
                                           upr = quantile(importance,0.975,na.rm=TRUE)),
                                        .(model,feature,seed)] %>%
        .[order(mean)] %>%
    mutate(significant = lwr>0 | upr<0)
})
names(mccv_importance_summary) <- dists

mccv_permuted_performance_summary <- lapply(dists,function(dist){
    mccv_permuted_data[[dist]]$performance[,.(lwr = quantile(value,0.025,na.rm=TRUE),
                                     mean = mean(value,na.rm=TRUE),
                                     upr = quantile(value,0.975,na.rm=TRUE)),
                                  .(model,metric,seed)] %>%
        .[order(mean)] %>%
        mutate(significant = !(lwr<0.5 & upr>0.5))
})
names(mccv_permuted_performance_summary) <- dists

mccv_permuted_prediction_summary <- lapply(dists,function(dist_){
    mccv_permuted_data[[dist_]]$patient_prediction[,.(
        class0_mean = mean(y_proba[y_pred==0]),
        class1_mean = mean(y_proba[y_pred==1]),
        fc = mean(y_proba[y_pred==1])/mean(y_proba[y_pred==0]),
        t_statistic = t.test(y_proba[y_pred==1],
                             y_proba[y_pred==0])[1],
        t_pvalue = t.test(y_proba[y_pred==1],
                          y_proba[y_pred==0])[3]
    ),
    .(model,seed)] %>%
        .[order(fc,decreasing = TRUE)] %>%
        mutate(significant = !(fc>2))
})
names(mccv_permuted_prediction_summary) <- dists

mccv_permuted_importance_summary <- lapply(dists,function(dist){
    mccv_permuted_data[[dist]]$feature_importance[,.(lwr = quantile(importance,0.025,na.rm=TRUE),
                                            mean = mean(importance,na.rm=TRUE),
                                            upr = quantile(importance,0.975,na.rm=TRUE)),
                                         .(model,feature,seed)] %>%
        .[order(mean)] %>%
        mutate(significant = !(lwr<0 & upr>0))
})
names(mccv_permuted_importance_summary) <- dists


# `perf_summary`, `pred_summary`, `imp_summary`, `permuted_perf_summary`, `permuted_pred_summary`, `permuted_imp_summary` --------------------------------------------------------------------

perf_summary <- function(seed_,digits_,metric_,model_){
    mccv_performance_summary[[dist_]][
        model==model_ &
            metric==metric_ &
            seed==seed_] %>%
        .[,paste0(
            signif(mean,digits_)," (",
            signif(lwr,digits_),", ",
            signif(upr,digits_),")")
        ]
}

pred_summary <- function(seed_,digits_,model_){
    mccv_prediction_summary[[dist_]][
        model==model_ &
            seed==seed_] %>%
        .[,paste0(
            signif(class1_mean,digits_)," vs. ",
            signif(class0_mean,digits_))
        ]
}

imp_summary <- function(seed_,digits_,model_){
    mccv_importance_summary[[dist_]][
        model==model_ &
            feature=="biomarker" &
            seed==seed_] %>%
        .[,paste0(
            signif(mean,digits_)," (",
            signif(lwr,digits_),", ",
            signif(upr,digits_),")")
        ]
}

permuted_perf_summary <- function(seed_,digits_,metric_,model_){
    mccv_permuted_performance_summary[[dist_]][
        model==model_ &
            metric==metric_ &
            seed==seed_] %>%
        .[,paste0(
            signif(mean,digits_)," (",
            signif(lwr,digits_),", ",
            signif(upr,digits_),")")
        ]
}

permuted_pred_summary <- function(seed_,digits_,model_){
    mccv_permuted_prediction_summary[[dist_]][
        model==model_ &
            seed==seed_] %>%
        .[,paste0(
            signif(class1_mean,digits_)," vs. ",
            signif(class0_mean,digits_))
        ]
}

permuted_imp_summary <- function(seed_,digits_,model_){
    mccv_permuted_importance_summary[[dist_]][
        model==model_ &
            feature=="biomarker" &
            seed==seed_] %>%
        .[,paste0(
            signif(mean,digits_)," (",
            signif(lwr,digits_),", ",
            signif(upr,digits_),")")
        ]
}

# Biomarker expression and prediction plots -----------------------------------------------------------

dir.create(paste0(img_dir,"biomarker_expression_prediction_plots"))
lapply(paste0(img_dir,"biomarker_expression_prediction_plots/",dists),dir.create) %>% invisible()
model_="Logistic Regression"
digits_=2
seeds=unique(datasets$normal$seed)
map2(rep(dists,500),rep(seeds,3),function(dist_,seed_){

    g1 <- mccv_biomarker_boxplot(seed_,dist = dist_) +
        ggtitle(
            paste0(
                str_to_title(dist_)," distributed data\n",
                "T-statistic: ",datasets_summary[[dist_]][seed==seed_,signif(t_statistic,2)],
                "\nP-value: ",datasets_summary[[dist_]][seed==seed_,signif(t_pvalue,2)]
            )
        )
    g2 <- mccv_prediction_boxplot(seed_,dist=dist_,model=model_) +
        ggtitle(
            paste0("AUROC: ",perf_summary(seed_,digits_,metric_ = "roc_auc",model_),"\n",
                   "Precision: ",perf_summary(seed_,digits_,metric_ = "average_precision",model_),"\n",
                   paste0("Log Odds: ",imp_summary(seed_,digits_,model_)))
        )

    cowplot::save_plot(
        paste0(
            img_dir,"biomarker_expression_prediction_plots/",
            dist_,"/",seed_,"_",
            janitor::make_clean_names(model_),
            ".pdf"),
        cowplot::plot_grid(g1,g2),
        base_width = 9,base_height=5
    )

}) %>% invisible()

# Significant biomarkers --------------------------------------------------

significant_performance_seeds <- lapply(dists,function(dist){
    mccv_performance_summary[[dist]][,
                                     .(N=sum(significant==TRUE)),
                                     .(model,seed)][
                                         N==mccv_performance_summary[[dist]][,n_distinct(metric)]
                                         ][,
        .(model,seed)] %>%
        unique()
})
names(significant_performance_seeds) <- dists

significant_prediction_seeds <- lapply(dists,function(dist){
    mccv_prediction_summary[[dist]][
        significant==TRUE,
        .(model,seed)] %>%
        unique()
})
names(significant_prediction_seeds) <- dists

significant_importance_seeds <- lapply(dists,function(dist){
    mccv_importance_summary[[dist]][
        significant==TRUE,
        .(model,seed)
        ] %>%
        unique()
})
names(significant_importance_seeds) <- dists

nonsignificant_performance_seeds <- lapply(dists,function(dist){
    mccv_performance_summary[[dist]][,
                                     .(N=sum(significant==FALSE)),
                                     .(model,seed)][
                                         N==mccv_performance_summary[[dist]][,n_distinct(metric)]
                                     ][,
                                       .(model,seed)] %>%
        unique()
})
names(nonsignificant_performance_seeds) <- dists

nonsignificant_prediction_seeds <- lapply(dists,function(dist){
    mccv_permuted_prediction_summary[[dist]][
        significant==FALSE,
        .(model,seed)] %>%
        unique()
})
names(nonsignificant_prediction_seeds) <- dists

nonsignificant_importance_seeds <- lapply(dists,function(dist){
    mccv_permuted_importance_summary[[dist]][
        significant==FALSE,
        .(model,seed)
        ] %>%
        unique() %>%
        .[order(model,seed)]
})
names(nonsignificant_importance_seeds) <- dists

significant_seeds <- lapply(dists,function(dist){
    merge(
        significant_performance_seeds[[dist]],
        significant_importance_seeds[[dist]],
        by=c("model","seed"),
        all=F
    ) %>%
        merge(
            significant_prediction_seeds[[dist]],
            by=c("model","seed"),
            all=F
        ) %>%
        unique() %>%
        .[order(model,seed)]
})
names(significant_seeds) <- dists

nonsignificant_seeds <- lapply(dists,function(dist){
    merge(
        nonsignificant_performance_seeds[[dist]],
        nonsignificant_importance_seeds[[dist]],
        by=c("model","seed"),
        all=T
    ) %>%
        merge(
            nonsignificant_prediction_seeds[[dist]],
            by=c("model","seed"),
            all=T
        ) %>%
        unique() %>%
        .[order(model,seed)]
})
names(nonsignificant_seeds) <- dists

really_significant_seeds <- lapply(dists,function(dist){
    merge(
        significant_seeds[[dist]],
        nonsignificant_seeds[[dist]],
        by=c("model","seed"),
        all=F
        ) %>%
        unique() %>%
        .[order(seed)]
})
names(really_significant_seeds) <- dists

map(really_significant_seeds,nrow)


# Significance venn diagrams (if more than one model) ----------------------------------------------

# for(dist_ in dists){
#     tmp <-
#         really_significant_seeds[[dist_]] %>%
#         mutate(mem=TRUE) %>%
#         pivot_wider(id_cols = "seed",
#                     names_from = "model",
#                     values_from = "mem",
#                     values_fill = FALSE)
#
#     g <- ggplot(tmp) +
#         ggvenn::geom_venn(aes(
#             A=`Logistic Regression`,
#             B=`Support Vector Machine`,
#             C=`Gradient Boosting Classifier`
#         ),
#         fill_color = colorspace::qualitative_hcl(3, palette = "Dark 3"),
#         text_size = 4
#         ) +
#         coord_fixed() +
#         theme_void(base_size = 16)
#
#     ggsave(paste0("docs/imgs/",dist_,"_distribution_model_significance_venn.pdf"),width=7,height=7)
#
# }

# Biomarker Volcano Plots ------------------------------------------------------------

lst <-
    lapply(dists,function(dist_){
        merge(
            mccv_performance_summary[[dist_]][,.(model,metric,seed,auroc = mean)],
            mccv_importance_summary[[dist_]][feature=="biomarker",.(model,seed,imp = mean)],
            by=c("model","seed"),
            all=T
        )
    })
names(lst) <- dists
g <-
    lst %>%
    bind_rows(.id="distribution") %>%
    filter((imp> -20 & imp<20)) %>%
    ggplot(aes(imp,auroc,fill=model)) +
    geom_point(pch=21,size=3) +
    xlab("Biomarker association") +
    ylab("Biomarker prediction (AUROC)") +
    facet_grid(model~metric) +
    scale_fill_brewer(palette = "Set1") +
    mytheme
ggsave(paste0(img_dir,"biomarker_association_vs_prediction.pdf"),width=6,height=4)

g <-
    lst %>%
    bind_rows(.id="distribution") %>%
    merge(
        datasets_summary %>% bind_rows(.id="distribution"),
        by=c("seed","distribution")
    ) %>%
    filter((fc> -20 & fc<20)) %>%
    ggplot(aes(fc,imp,fill=auroc)) +
    geom_point(pch=21,size=3) +
    xlab("Biomarker Fold Change") +
    ylab("Biomarker Association") +
    guides(fill=guide_colorbar(title="Biomarker prediction (AUROC)",
                               title.position = "top",
                               barwidth = 15)) +
    facet_wrap(~model,nrow=2) +
    scale_fill_gradientn(colours = c("cyan", "white", "red"),
                         values = scales::rescale(c(0.3, 0.8,1))) +
    ggthemes::theme_clean(base_size = 16) +
    theme(
        legend.position = "bottom"
    )
ggsave(paste0(img_dir,"biomarker_foldchange_vs_association.pdf"),width=5,height=4)

g <-
    lst %>%
    bind_rows(.id="distribution") %>%
    merge(
        datasets_summary %>% bind_rows(.id="distribution"),
        by=c("seed","distribution")
    ) %>%
    filter((fc> -20 & fc<20)) %>%
    ggplot(aes(fc,auroc,fill=model)) +
    geom_point(pch=21,size=3) +
    geom_vline(xintercept = 1,color="black",size=2) +
    xlab("Biomarker Fold Change") +
    ylab("Biomarker Prediction (AUROC)") +
    scale_fill_brewer(palette = "Set1") +
    guides(fill=guide_legend(title.position = "top")) +
    facet_wrap(~model,nrow=2) +
    mytheme +
    labs(caption="x-intercept = 1")
ggsave(paste0(img_dir,"biomarker_foldchange_vs_prediction.pdf"),width=6,height=4)

g <-
    lst %>%
    bind_rows(.id="distribution") %>%
    merge(
        datasets_summary %>% bind_rows(.id="distribution"),
        by=c("seed","distribution")
    ) %>%
    filter((fc> -20 & fc<20)) %>%
    ggplot(aes(log2fc,auroc,fill=model)) +
    geom_point(pch=21,size=3) +
    xlab("Biomarker Log2FoldChange") +
    ylab("Biomarker Prediction (AUROC)") +
    scale_fill_brewer(palette = "Set1") +
    guides(fill=guide_legend(title.position = "top")) +
    facet_wrap(~model,nrow=2) +
    mytheme
ggsave(paste0(img_dir,"biomarker_log2foldchange_vs_prediction.pdf"),width=6,height=4)

g <-
    lst %>%
    bind_rows(.id="distribution") %>%
    merge(
        datasets_summary %>% bind_rows(.id="distribution"),
        by=c("seed","distribution")
    ) %>%
    mutate(t_padj = p.adjust(t_pvalue,"fdr")) %>%
    merge(
        really_significant_seeds %>%
            bind_rows(.id="distribution") %>%
            mutate(significant=TRUE),
        by=c("seed","distribution","model"),
        all.x=T
    ) %>%
    mutate(
        significant = ifelse(is.na(significant),FALSE,TRUE)
    ) %>%
    filter((fc> -20 & fc<20)) %>%
    ggplot(aes(-log10(t_pvalue),auroc,fill=significant)) +
    geom_point(pch=21,size=3) +
    scale_x_sqrt(labels=c("0.05","FWER","10","20","40","60","80"),
                 breaks=-log10(
                     c(0.05,1/1500,1E-10,1E-20,1E-40,1E-60,1E-80)
                 )) +
    geom_vline(xintercept = -log10(0.05),color="black",size=2) +
    geom_vline(xintercept = -log10(1/1500),color="black",size=2) +
    xlab("Biomarker t-test P-value") +
    ylab("Biomarker Prediction (AUROC)") +
    scale_fill_brewer(palette = "Set1",direction = -1) +
    guides(fill=guide_legend(title="Predictive?",title.position = "top")) +
    facet_wrap(~model,nrow=2) +
    mytheme +
    theme(
        axis.text.x = element_text(angle=45,vjust=1,hjust=1),
        legend.position = "bottom"
    )
ggsave(paste0(img_dir,"biomarker_pvalue_vs_prediction_by_significance.pdf"),width=8,height=6)


# Normal distribution data performance ------------------------------------

model_="Logistic Regression"
map(c("roc_auc","average_precision"),function(metric_){
    g <- lapply(dists,function(dist_){
    datasets[[dist_]] %>%
        dcast(seed + distribution ~ class,
              value.var = "biomarker",
              fun.aggregate = mean)}) %>%
    bind_rows() %>%
    rename(class0=`0`,class1=`1`) %>%
    .[order(seed,distribution)] %>%
    merge(
        mccv_performance_summary %>%
            bind_rows(.id="distribution") %>%
            .[model==model_ & metric==metric_,
              .(seed,distribution,auroc=mean)] %>%
            .[order(seed,distribution)],
        all=T,
        by=c("seed","distribution")
    ) %>%
    ggplot(aes(class0,class1,fill=auroc)) +
    geom_point(pch=21,size=3,color="black") +
    geom_abline(slope=1,intercept=0,linetype="dashed",color="red") +
    colorspace::scale_fill_continuous_diverging(
        palette = "Blue-Red",mid=mccv_performance_summary[[dist_]][,mean(mean)]) +
    guides(fill=guide_colorbar(title=metric_,title.position = "top",
                               barwidth = 15)) +
    ggthemes::theme_clean(base_size = 16) +
    xlab("Class 0 biomarker expression") +
    ylab("Class 1 biomarker expression") +
    theme(
        legend.position = "bottom"
    )
ggsave(paste0("docs/imgs/",metric_,"_prediction_by_average_biomarker_class_difference.pdf"),width=5,height=5)
})
