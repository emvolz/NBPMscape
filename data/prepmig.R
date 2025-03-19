# Format migration data for inclusion in jul code 


K = readRDS('./ITL2_key.rds' ) 
M = readRDS('./commuting_ITL2_matrix.rds') 
# K1 <- K[ grepl( K$code, patt= '^TL.*' ), ] 

M1 <- M[ -c(2,3), ]
rownames(M1)[1] <- 'na'
cnms = colnames(M1) 
M1 <- M1[ c('na', cnms), cnms]
m1l = lapply( 1:ncol(M1), function(k) M1[,k] )
names(m1l) <- cnms
p1l <- lapply( m1l, function(x) setNames( x / sum(x) , names(x) )) |> setNames( names(m1l) )

b1l <- lapply( cnms, function(nm) M1[nm,]/sum(M1[nm,]) ) |> setNames( cnms)

K1 <- K[ !duplicated(K$code), ]
rownames(K1) <- K1$code  
K1 <- K1[ cnms, ]

saveRDS(K1, './ITL2_key2.rds') 
saveRDS( m1l, './commuting_ITL2_list.rds')
saveRDS( p1l, './commuting_ITL2_prob_list.rds')
saveRDS( b1l, './commuting_ITL2_inprob_list.rds')

# p1l[[1]] |> plot()
