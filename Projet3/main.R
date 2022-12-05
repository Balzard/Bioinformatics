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
i
