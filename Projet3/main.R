train = read.csv("./Projet3/train.csv.bz2", header = TRUE)
train = data.frame(train[,-1], row.names = train[,1]) # to get patients IDs as row name

ns_filtering = function(df, kept){
    df = subset(df, select = -labels)
    vec = c(apply(df, 2, var))
    names(vec) = names(df)
    vec = unlist(vec)
    vec = sort(vec, decreasing = TRUE)
    vec[1:(length(vec)*kept)]
}

order_train = ns_filtering(train, 0.25)
length(order_train)
names(order_train)
saveRDS(order_train, file = "./Projet3/probesets_largest_var.rds")


subset_df = function(df){
    filter_features = names(ns_filtering(df, 0.25))
    filter_features = append(filter_features, "labels")
    df = subset(df, select = filter_features)
}

differential_selection = function(data, target, alpha){
    features = names(data)[1:(length(names(data)) - 1)]
    counter = 0
    for(feature in features){
        feature1 = subset(data[data[, target] == "No tumor",], select = feature)
        feature2 = subset(data[data[, target] == "Glioblastoma",], select = feature)
        test = t.test(feature1, feature2)
        p_value = test$p.value
        if(p_value <= alpha){
            counter = counter + 1
        }
    }
    counter
}

subset_train = subset_df(train)
i = differential_selection(subset_train, "labels", 0.05)


correct = function(data, target, alpha, method){
    features = names(data)[1:(length(names(data)) - 1)]
    counter = 0
    nb_tests = length(features)
    bonferroni_correction = alpha / nb_tests
    p_values = vector()
    bonferroni_features = vector()

    for(feature in features){
        feature1 = subset(data[data[, target] == "No tumor",], select = feature)
        feature2 = subset(data[data[, target] == "Glioblastoma",], select = feature)
        test = t.test(feature1, feature2)
        p_value = test$p.value

        if(method == "bonferroni" && p_value <= bonferroni_correction){
            p_values = append(p_values, p_value)
            bonferroni_features = append(bonferroni_features, feature)
        }

        if(method == "fdr"){
            p_values = append(p_values, p_value)
        }
    }
    if(method == "bonferroni"){
        names(p_values) = bonferroni_features
        p_values = sort(p_values)
        return(p_values)
    }
    if(method == "fdr"){
        names(p_values) = features
        p_values = sort(p_values)
        for(i in 1:length(p_values)){
            test = (p_values[i] * nb_tests) / i
            if(test >= alpha){
                p_values = p_values[1:(i - 1)]
                return(p_values)
            }
        }
    }
}

bonferroni_values = correct(subset_train, "labels", 0.05, "bonferroni")
saveRDS(bonferroni_values, "./Projet3/bonferroni_values.rds")
fdr_values = correct(subset_train, "labels", 0.05, "fdr")
saveRDS(fdr_values, "./Projet3/fdr_values.rds")


getMostDiffFeatures = function(data, target, limit=50){
    features = names(data)[1:(length(names(data)) - 1)]
    p_values = vector()
    for(feature in features){
        feature1 = subset(data[data[, target] == "No tumor",], select = feature)
        feature2 = subset(data[data[, target] == "Glioblastoma",], select = feature)
        test = t.test(feature1, feature2)
        p_value = test$p.value
        p_values = append(p_values, p_value)
    }
    names(p_values) = features
    p_values = sort(p_values)
    return(p_values[1:limit])
}

most_diff_features = getMostDiffFeatures(subset_train, "labels")
saveRDS(most_diff_features, "./Projet3/most_diff_features.rds")


library("gplots")

getHeatmap = function(subset_train, most_diff_features){
    train_most_diff_features = subset(subset_train, select = names(most_diff_features))
    matrix_heatmap = data.matrix(train_most_diff_features)
    labels_color = c("No tumor"="blue", "Glioblastoma"="purple")
    PID_colors = vector()

    for(i in train[, "labels"]){
        PID_colors = append(PID_colors, labels_color[i])
    }

    print(length(matrix_heatmap[,1:50]))

    png(file = "./Projet3/heatmap.png", width = 800, height = 1000)
    heatmap.2(x=matrix_heatmap, trace = "none", main = "Heatmap of features and samples", xlab = "Genes", ylab = "Patients ID", RowSideColors = PID_colors)
    legend(x="topright",legend = c("No tumor", "Glioblastoma"), fill = c("blue", "purple"))
    dev.off()

}

most_diff_features = readRDS("./Projet3/most_diff_features.rds")
most_diff_features

getHeatmap(subset_train, most_diff_features)

most_diff_features
