library(HydeNet)
#make a df
mtcars2 <- transform(mtcars,
                     cyl = factor(cyl),
                     gear=factor(gear),
                     am = factor(am))

#Use the HydeNetwork() function to set up the network structure. Passed the training data frame mtcars2 
#to be used for populating the parameters of the network.
#The method for HydeNetwork() uses syntax similar to the dag() command in the gRbase library.
carNet <- HydeNetwork(~ cyl
                      + disp | cyl
                      + hp | disp
                      + wt
                      + gear
                      + mpg | disp*hp*wt*gear,
                      data=mtcars2)

# the above is the automatic way to calc cpd of net. cpds of categorical nodes with all 
#categorical parents are estimated using tabulation (see help('cpt')), while conditional 
#distributions for categorical nodes with at least one continuous parent node are 
#estimated using glm(, family="binomial").

#can also manaully set it up using
#by calling HydeNetwork() with a list of model objects or specifying only the network structure 
#in a call to HydeNetwork() and then calling setNode() and setNodeModels() to 
#specify distributions for each individual node.

#see “Working with HydeNet Objects” vignette (vignette("WorkingWithHydeNetObjects", package = "HydeNet"))
#for the details of HydeNetwork(), setNode() and setNodeModels().

#Use plot method for HydeNetwork objects.

plot(carNet)

#It using the DiagrammeR library. There are different default node shapes and colors for 
#random variable nodes, deterministic nodes, decision nodes and utility nodes. 
#See  “Building and Customizing HydeNet Plots” vignette (vignette("HydeNetPlots", package = "HydeNet"))
#for details on customization.

#HydeNet interfaces with JAGS through the rjags package. Given a fully-specified 
#HydeNetwork model, the function writeNetworkModel() will create the JAGS model 
#script (by calling the writeJagsModel() function for each node, which creates 
#each node’s portion of the script). Here’s an example (note, as indicated by the need to 
#use a triple colon, that writeJagsModel() is not an exported function to the HydeNet 
#namespace):

HydeNet:::writeJagsModel(carNet, node = "cyl")
HydeNet:::writeJagsModel(carNet, node = "mpg")
writeNetworkModel(carNet, pretty = TRUE)

#The function compileJagsModel() (or compileDecisionModel(), if the network 
#contains decision and utility nodes) is called. This function is a wrapper for the 
#jags.model() function in the rjags library which, in addition to the JAGS model script, 
#passes any conditional probability tables from the network model into JAGS as arrays. 
#(Note: all other arguments to jags.model(), such as the number of chains initial values, 
#and the number of iterations for adaptation are passed through to the function.)

#The functions compileJagsModel() and compileDecisionModel() essentially return augmented 
#objects of class jags. We label these objects with the class compiledHydeNetwork.

carNet1 <- compileJagsModel(carNet)

#Evidence on individual nodes can also be entered with compileJagsModel(), thus 
#facilitating dynamic inference as data are observed. See the description of the data 
#argument in help(rjags::jags.model) for details.

carNet2 <- compileJagsModel(carNet, data = list(cyl = "8") )

#As a convenience, HydeNet can also translate labels into the numeric codes

carNet3 <- compileJagsModel(carNet, data=list(cyl="8"))

#MCMC sampling is performed by calling HydePosterior(). This function calls 
#coda.samples() from the rjags library and returns an augmented mcmc.list object. 
#MCMC diagnostics and inference from this point onward are performed in the same 
#manner as would typically be done when using rjags.

post1 <- HydePosterior(carNet1,
                       variable.names = c("cyl","hp","mpg"),
                       n.iter = 10000,
                       bind=FALSE)

post2 <- HydePosterior(carNet2,
                       variable.names = c("cyl","hp","mpg"),
                       n.iter = 10000,
                       bind = FALSE)
str(post1, max.level = 3)

plot(post1$codas[,c("hp","mpg")])
