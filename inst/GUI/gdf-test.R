w <- gwindow("gdf test")
g <- ggroup(horizontal=FALSE, cont=w)
d <- gdf(data.frame(a=c(1.1, 2.2, 3.3), b= c("one", "two", "three"), stringsAsFactors=FALSE), 
  name = "Test for editing character columns", cont=g)
size(d) <- c(400,400)
gbutton("Change [1,2]", handler=function(h,...) {
  d[1,2] <- "test"
}, container=g)
gbutton("Value of [1,2]", handler=function(h,...) {
  galert(d[1,2], parent=d)
}, cont=g)
visible(w) <- TRUE
