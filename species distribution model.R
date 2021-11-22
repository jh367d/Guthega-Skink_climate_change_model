

r <- getOption("repos")
r["CRAN"] <- "http://cran.ms.unimelb.edu.au/"
options(repos=r)
# print warnings immediately
options(warn=1)

# print out list of installed packages
write.table(installed.packages()[,c("Package", "Version", "Priority")],
            row.names=FALSE)

#script to run to develop distribution models
###check if libraries are installed, install if necessary and then load them
necessary=c("ggplot2","tools", "rjson", "dismo","SDMTools", "gbm", "rgdal", "rgeos", "pROC", "png", "gstat", "biomod2", "gdalUtils", "spatial.tools") #list the libraries needed
installed = necessary %in% installed.packages() #check if library is installed
if (length(necessary[!installed]) >=1) {
    install.packages(necessary[!installed], dep = T) #if library is not installed, install it
}
for (lib in necessary) {
    library(lib,character.only=T) #load the libraries
}

# load parameters
params = rjson::fromJSON(file="params.json")
bccvl.params <- params$params
bccvl.env <- params$env
rm(params)

# set working directory (script runner takes care of it)
setwd(bccvl.env$outputdir)
# Set raster tmpdir - we do this because raster sometimes makes
# temp files (e.g. when cropping).
# Might want to make this configurable - e.g. we might want to
# control maxmemory, and/or other raster options
rasterOptions(tmpdir=paste(bccvl.env$workdir,"raster_tmp",sep="/"))

# Use seed supplied if any. Otherwise generate a random seed.
seed = bccvl.params$random_seed
if (is.null(seed)) {
    seed = runif(1, -2^31, 2^31-1)
}
seed = as.integer(seed)
set.seed(seed)
bccvl.params["random_seed"] = seed


############################################################
#
# define helper functions to print parameters
#
############################################################

parameter.as.string <- function (param, value) {
    pname <- gsub("_", " ", param)
    if (param == "prevalence") {
        pname = "weighted response weights"
    }
    else if (param == "var_import") {
        pname = "resampling"
    }
    else if (param == "nbcv") {
        pname = "NbCV"
    }
    else if (param == "n_trees") {
        pname = "trees added each cycle"
    }
    else if (param == "control_xval") {
        pname = "cross-validations"
    }
    else if (param == "control_minbucket") {
        pname = "minimum bucket"
    }
    else if (param == "control_minsplit") {
        pname = "minimum split"
    }
    else if (param == "control_cp") {
        pname = "complexity parameter"
    }
    else if (param == "control_maxdepth") {
        pname = "maximum depth"
    }
    else if (param == "irls_reg") {
        pname = "irls.reg"
    }
    else if (param == "maxit") {
        pname = "maximum iterations"
    }
    else if (param == "mgcv_tol") {
        pname = "convergence tolerance"
    }
    else if (param == "mgcv_half") {
        pname = "number of halvings"
    }
    else if (param == "n_minobsinnode") {
        pname = "Min observations in terminal node"
    }
    else if (param == "control_epsilon") {
        pname = "control: epsilon"
    }
    else if (param == "control_maxit") {
        pname = "control: maxit"
    }
    else if (param == "control_trace") {
        pname = "control: trace"
    }
    else if (param == "model") {
        pname = "Model returned"
    }
    else if (param == "x") {
        pname = "x returned"
    }
    else if (param == "y") {
        pname = "y returned"
    }
    else if (param == "qr") {
        pname = "QR returned"
    }
    else if (param == "singular_ok") {
        pname = "Singular fit ok"
    }
    else if (param == "thresh") {
        pname = "threshold"
    }
    else if (param == "maximumiterations") {
        pname = "Maximum iterations"
    }
    else if (param == "ntree") {
        pname = "number of trees"
    }
    else if (param == "mtry") {
        pname = "number of variables at each split (mtry)"
    }
    else if (param == "nodesize") {
        pname = "node size"
    }
    else if (param == "maxnodes") {
        pname = "maximum nodes"
    }
    else if (param == "pa_ratio") {
        pname = "absence-presence ratio"
    }
    return(paste(pname, " = ", value, "\n", sep="", collapse=""))
}

parameter.print <- function(params) {
    func = params[["function"]]
    if (is.null(func))
        return("")
    cat("Algorithm:", func, "\n")

    pnames = c("random_seed")
    if (func == "ann") {
        pnames = c("prevalence", "var_import", "maxit", "nbcv", "rang", "random_seed")
    }
    else if (func == "brt") {
        pnames = c("tree_complexity", "learning_rate", "bag_fraction", "n_folds", "prev_stratify", "family", "n_trees", "max_trees", "tolerance_method", "tolerance_value", "random_seed")
    }
    else if (func == "cta") {
        pnames = c("prevalence", "var_import", "method", "control_xval", "control_minbucket", "control_minsplit", "control_cp", "control_maxdepth", "random_seed")
    }
    else if (func == "fda") {
        pnames = c("prevalence", "var_import", "method", "random_seed")
    }
    else if (func == "gam") {
        pnames = c("prevalence", "var_import", "interaction_level", "family", "irls_reg", "epsilon", "maxit", "mgcv_tol", "mgcv_half", "random_seed")
    }
    else if (func == "gamlss") {
        pnames = c("sigma_formula", "nu_formula", "tau_formula", "family", "weights", "contrasts", "method", "start_from", "mu_start", "sigma_start", "nu_start", "tau_start", "mu_fix", "sigma_fix", "nu_fix", "tau_fix", "control", "i_control", "other_args", "random_seed")
    }
    else if (func == "gbm") {
        pnames = c("prevalence", "var_import", "distribution", "n_trees", "interaction_depth", "n_minobsinnode", "shrinkage", "bag_fraction", "train_fraction", "cv_folds", "random_seed")
    }
    else if (func == "glm") {
        pnames = c("prevalence", "var_import", "type", "interaction_level", "test", "family", "mustart", "control_epsilon", "control_maxit", "control_trace", "random_seed")
    }
    else if (func == "lm") {
        pnames = c("subset", "weights", "na_action", "method", "model", "x", "y", "qr", "singular_ok", "contrasts", "offset", "random_seed")
    }
    else if (func == "manova") {
        pnames = c("projections_returned", "qr", "contrasts", "subset", "weights", "na_action", "random_seed")
    }
    else if (func == "mars") {
        pnames = c("prevalence", "var_import", "degree", "nk", "penalty", "thresh", "prune", "random_seed")
    }
    else if (func == "maxent") {
        pnames = c("prevalence", "var_import", "maximumiterations", "linear", "quadratic", "product", "threshold", "hinge", "lq2lqptthreshold", "lq2lqthreshold", "hingethreshold", "beta_threshold", "beta_categorical", "beta_lqp", "beta_hinge", "defaultprevalence", "random_seed")
    }
    else if (func == "rf") {
        pnames = c("prevalence", "var_import", "do.classif", "ntree", "mtry", "nodesize", "maxnodes", "random_seed")
    }
    else if (func == "sre") {
        pnames = c("prevalence", "var_import", "quant", "random_seed")
    }

    pnames = c(pnames, "pa_ratio", "pa_strategy", "pa_sre_quant", "pa_disk_min", "pa_disk_max")
    for (p in pnames) {
        cat(parameter.as.string(p, params[[p]]))
    }
    return("")
}


# Print out parameters used
parameter.print(bccvl.params)


############################################################
#
# define helper functions to use in bccvl
#
############################################################

## Needed for tryCatch'ing:
bccvl.err.null <- function (e) return(NULL)

# read species presence/absence data
#    return NULL if filename is not  given
# TODO: shall we set projection here as well? use SpatialPoints?
bccvl.species.read <- function(filename, month_filter=NULL) {
    if (!is.null(filename)) {
        # We might loose precision of lon/lat when ronverting to double,
        # However, given the nature of the numbers, and the resolution of raster files
        # we deal with, this shouldn't be a problem.
        csvfile = read.csv(filename, colClasses=c("lon"="numeric", "lat"="numeric"))

        # keep only lon and lat columns; for MM include month column
        if (is.null(month_filter)) {
            csvfile = csvfile[c("lon","lat")]
            return(csvfile)
        }
        else {
            csvfile = csvfile[c("lon","lat","month")]
            return(subset(csvfile, month %in% unlist(month_filter)))
        }
    }
}

bccvl.data.transform <- function(data, climate.data)
{
    if (!is.null(data) & !compareCRS(data, climate.data, verbatim=TRUE)) {
        sp <- SpatialPoints(data)
        if (is.na(crs(sp))) {
            crs(sp) <- '+init=epsg:4326'
        }

        newdata <- as.data.frame(spTransform(sp, crs(climate.data)))
        names(newdata) <- names(data)
        return(newdata)
    }
    return(data)
}

bccvl.format.outfilename <- function(filename, id_str, ext)
{
    return(sprintf("%s_%s.%s", filename, id_str, ext))
}

# BIOMOD_FormatingData(resp.var, expl.var, resp.xy = NULL, resp.name = NULL, eval.resp.var = NULL,
#   eval.expl.var = NULL, eval.resp.xy = NULL, PA.nb.rep = 0, PA.nb.absences = 1000, PA.strategy = 'random',
#   PA.dist.min = 0, PA.dist.max = NULL, PA.sre.quant = 0.025, PA.table = NULL, na.rm = TRUE)
#
# resp.var a vector, SpatialPointsDataFrame (or SpatialPoints if you work with `only presences' data) containing species data (a single species) in binary format (ones for presences, zeros for true absences and NA for indeterminated ) that will be used to build the species distribution models.
# expl.var a matrix, data.frame, SpatialPointsDataFrame or RasterStack containing your explanatory variables that will be used to build your models.
# resp.xy optional 2 columns matrix containing the X and Y coordinates of resp.var (only consider if resp.var is a vector) that will be used to build your models.
# eval.resp.var a vector, SpatialPointsDataFrame your species data (a single species) in binary format (ones for presences, zeros for true absences and NA for indeterminated ) that will be used to evaluate the models with independant data (or past data for instance).
# eval.expl.var a matrix, data.frame, SpatialPointsDataFrame or RasterStack containing your explanatory variables that will be used to evaluate the models with independant data (or past data for instance).
# eval.resp.xy opional 2 columns matrix containing the X and Y coordinates of resp.var (only consider if resp.var is a vector) that will be used to evaluate the modelswith independant data (or past data for instance).
# resp.name response variable name (character). The species name.
# PA.nb.rep number of required Pseudo Absences selection (if needed). 0 by Default.
# PA.nb.absences number of pseudo-absence selected for each repetition (when PA.nb.rep > 0) of the selection (true absences included)
# PA.strategy strategy for selecting the Pseudo Absences (must be `random', `sre', `disk' or `user.defined')
# PA.dist.min minimal distance to presences for `disk' Pseudo Absences selection (in meters if the explanatory is a not projected raster (+proj=longlat) and in map units (typically also meters) when it is projected or when explanatory variables are stored within table )
# PA.dist.max maximal distance to presences for `disk' Pseudo Absences selection(in meters if the explanatory is a not projected raster (+proj=longlat) and in map units (typically also meters) when it is projected or when explanatory variables are stored within table )
# PA.sre.quant quantile used for `sre' Pseudo Absences selection
# PA.table a matrix (or a data.frame) having as many rows than resp.var values. Each column correspund to a Pseudo-absences selection. It contains TRUE or FALSE indicating which values of resp.var will be considered to build models. It must be used with `user.defined' PA.strategy.
# na.rm logical, if TRUE, all points having one or several missing value for environmental data will be removed from analysis


# This uses the biomods function BIOMOD_FormatingData to format user input data.
# It generates pseudo absence points if true absence data are not available or
# adds pseudo absence data to an existing absence dataset.
bccvl.biomod2.formatData <- function(true.absen=NULL,
                                  pseudo.absen.points=0,
                                  pseudo.absen.strategy='random',
                                  pseudo.absen.disk.min=0,
                                  pseudo.absen.disk.max=NULL,
                                  pseudo.absen.sre.quant = 0.025,
                                  climate.data=NULL,
                                  occur=NULL,
                                  species.name=NULL,
                                  save.pseudo.absen=TRUE,
                                  save.env.absen=TRUE,
                                  save.env.occur=TRUE,
                                  generate.background.data=FALSE,
                                  species_algo_str=NULL) {

    # Initialise parameters to default value if not specified
    if (is.null(pseudo.absen.strategy)) {
        pseudo.absen.strategy = 'random'
    }
    if (is.null(pseudo.absen.disk.min)) {
        pseudo.absen.disk.min = 0
    }
    if (is.null(pseudo.absen.sre.quant)) {
        pseudo.absen.sre.quant = 0.025
    }

    # Read true absence point if available.
    if (is.null(true.absen)) {
        # create an empty data frame for bkgd points
        absen = data.frame(lon=numeric(0), lat=numeric(0))
        # To generate pseudo=absence points
        pseudo.absen.rep = 1
        if (!save.pseudo.absen) {
            pseudo.absen.rep = 0
        }
    }
    else {
        # Ensure true absence dataset is in same projection system as climate.
        absen <- true.absen
        if (!is.null(climate.data) && nrow(true.absen) > 0) {
            absen <- bccvl.data.transform(true.absen, climate.data)
        }

        # Do not generate pseudo absence point when true absence points are available
        pseudo.absen.rep = 0
        pseudo.absen.strategy = 'none'
        pseudo.absen.points = nrow(absen)
        cat("No pseudo absence point is generated.")
    }

    # Generate background data as pseudo absence points
    if (pseudo.absen.strategy != 'none' & generate.background.data) {
        biomod.data.pa <- c(rep(1, nrow(occur)), rep(0, nrow(absen)))
        myBackgrdData <-
            BIOMOD_FormatingData(resp.var  = biomod.data.pa,
                                 expl.var  = climate.data,
                                 resp.name = species.name,
                                 PA.nb.rep = pseudo.absen.rep,
                                 PA.nb.absences = pseudo.absen.points,
                                 PA.strategy = pseudo.absen.strategy,
                                 PA.dist.min = pseudo.absen.disk.min,
                                 PA.dist.max = pseudo.absen.disk.max,
                                 PA.sre.quant = pseudo.absen.sre.quant)

        # Get background data as absence data
        colnames(myBackgrdData@coord) <- c('lon', 'lat')
        absen <- myBackgrdData@coord[c(which(is.na(myBackgrdData@data.species))), c('lon', 'lat')]

        # Do not generate pa in next call to BIOMOD_FormatingData
        pseudo.absen.rep = 0
        pseudo.absen.strategy = 'none'
    }

    biomod.data <- rbind(occur[,c("lon", "lat")], absen[,c("lon", "lat")])
    biomod.data.pa <- c(rep(1, nrow(occur)), rep(0, nrow(absen)))

    myBiomodData <-
        BIOMOD_FormatingData(resp.var  = biomod.data.pa,
                             expl.var  = climate.data,
                             resp.xy   = biomod.data,
                             resp.name = species.name,
                             PA.nb.rep = pseudo.absen.rep,
                             PA.nb.absences = pseudo.absen.points,
                             PA.strategy = pseudo.absen.strategy,
                             PA.dist.min = pseudo.absen.disk.min,
                             PA.dist.max = pseudo.absen.disk.max,
                             PA.sre.quant = pseudo.absen.sre.quant)

    # Save the pseudo absence points generated to file
    pa_filename = bccvl.format.outfilename(filename="pseudo_absences", id_str=species_algo_str, ext="csv")
    absenv_filename = bccvl.format.outfilename(filename="absence_environmental", id_str=species_algo_str, ext="csv")
    occenv_filename = bccvl.format.outfilename(filename="occurrence_environmental", id_str=species_algo_str, ext="csv")
    if (pseudo.absen.rep != 0) {
        pseudoAbsen = myBiomodData@coord[c(which(is.na(myBiomodData@data.species))), c('lon', 'lat')]
        if (save.pseudo.absen & nrow(pseudoAbsen) > 0) {
            bccvl.write.csv(pseudoAbsen, pa_filename, rownames = FALSE)
        }

        # save the pseudo absence points with environmental variables
        if (save.env.absen) {
            bccvl.merge.save(climate.data, pseudoAbsen, species.name, absenv_filename)
        }
    }
    else if (nrow(absen) > 0) {
        # save true-absence/background data generated
        if (!is.null(true.absen)) {
            # rename true-absence file
            pa_filename = bccvl.format.outfilename(filename="absence", id_str=species_algo_str, ext="csv")
        }
        bccvl.write.csv(absen, pa_filename, rownames = FALSE)

        # save the true absence points/background points with environmental variables
        if (save.env.absen) {
            bccvl.merge.save(climate.data, absen, species.name, absenv_filename)
        }
    }

    # save the occurrence datasets with environmental variables
    if (save.env.occur) {
        bccvl.merge.save(climate.data, occur, species.name, occenv_filename)
    }

    return(myBiomodData)
}

bccvl.merge.save <- function(env, csvdata, spname, ofname)
{
  data = cbind(csvdata, species=spname, extract(env, csvdata[c('lon','lat')]))

  bccvl.write.csv(data, ofname, rownames=FALSE)
}

# warning was doing odd things. I just want to print the deng thing.
bccvl.log.warning <-function(str, prefix="BCCVL Warning: ")
{
    print(paste(prefix, str, sep=""))
}

bccvl.raster.load <- function(filename) {
    # load raster and assign crs if missing
    r = raster(filename)
    if (is.na(crs(r))) {
        crs(r) = CRS("+init=epsg:4326")
    }
    return(r)
}


# rasters: a vector of rasters, all rasters should have same resolution
# common.crs: crs to use to calculate intersection
bccvl.raster.common.extent <- function(rasters, common.crs)
{
    # bring all rasters into common crs
    extent.list = lapply(rasters, function(r) { extent(projectExtent(r, common.crs)) })
    # intersect all extents
    common.extent = Reduce(intersect, extent.list)
    # compare all against commen extents to find out if all extents are the same (used to print warning)
    equal.extents = all(sapply(extent.list, function (x) common.extent == x))

    return (list(equal.extents=equal.extents, common.extent=common.extent))
}

bccvl.raster.extent.to.str <- function(ext)
{
  return(sprintf("xmin=%f xmax=%f ymin=%f ymax=%f", ext@xmin, ext@xmax, ext@ymin, ext@ymax));
}

# rasters: a vector of rasters ... preferrably empty
# resamplingflag: a flag to determine which resampling approach to take
# selected_layers: a list of indexes to the raster layers to be considered when determine the resolution to be used.
bccvl.rasters.common.resolution <- function(rasters, resamplingflag, selected_layers) {
    # Return the resolution of the raster given by the index
    get.resolution <- function(i, rasters)
    {
        return(res(rasters[[i]]))
    }

    # Get resolutions of the raster layers
    if (is.null(selected_layers)) {
        resolutions = lapply(rasters, res)
    }
    else {
        # get the resolutions of the of the selected raster layers only 
        resolutions = lapply(selected_layers, get.resolution, rasters)
    }
    
    #resolutions = lapply(rasters, res)
    if (resamplingflag == "highest") {
        common.res = Reduce(pmin, resolutions)
    } else if (resamplingflag == "lowest") {
        common.res = Reduce(pmax, resolutions)
    }

    # Get resolutions of all input raster layers
    resolutions = lapply(rasters, res)
    is.same.res = all(sapply(resolutions, function(x) all(common.res == x)))
    return (list(common.res=common.res, is.same.res=is.same.res))
}

# generate reference raster with common resolutin, crs and extent
bccvl.rasters.common.reference <- function(rasters, resamplingflag, selected_layers) {
    # create list of empty rasters to speed up alignment
    empty.rasters = lapply(rasters, function(x) { projectExtent(x, crs(x)) })
    # choose a common.crs if all crs in rasters are the same use that one, otherwise use EPSG:4326 (common data in bccvl)
    common.crs = crs(empty.rasters[[1]])
    # TODO: print warning about reprojecting if necessary (if inside next condition)
    if (! do.call(compareRaster, c(empty.rasters, extent=FALSE, rowcol=FALSE, prj=TRUE, res=FALSE, orig=FALSE, rotation=FALSE, stopiffalse=FALSE))) {
        # we have different CRSs, so use EPSG:4326 as common
        # TODO: another strategy to find common CRS?
        common.crs = CRS("+init=epsg:4326")
        # project all rasters into common crs
        bccvl.log.warning(sprintf("Auto projecting to common CRS %s", common.crs))
        empty.rasters = lapply(empty.rasters, function(x) { projectExtent(x, common.crs) })
    }

    # determine commen.extent in common.crs
    # Note: extent is in projection units, -> rasters have to be in same CRS
    ce = bccvl.raster.common.extent(empty.rasters, common.crs)
    if (! ce$equal.extents) {
        bccvl.log.warning(sprintf("Auto cropping to common extent %s", bccvl.raster.extent.to.str(ce$common.extent)))
    }

    # determine common resolution
    # Note: resolution is usually in projection units. -> rasters should be in same CRS
    cr = bccvl.rasters.common.resolution(empty.rasters, resamplingflag, selected_layers)
    # TODO: print warning about resampling: common.res$is.same.res
    if (! cr$is.same.res) {
        bccvl.log.warning(sprintf("Auto resampling to %s resolution [%f %f]", resamplingflag, cr$common.res[[1]], cr$common.res[[2]]))
    }

    # apply common extent and resolution to empty rasters
    empty.rasters = lapply(
        empty.rasters,
        function(x) {
            extent(x) = ce$common.extent
            res(x) = cr$common.res
            return(x)
        })
    # from now an all empty.rasters should be exactly the same, pick first and return as
    # template.
    return(empty.rasters[[1]])
}

bccvl.rasters.warp <- function(raster.filenames, raster.types, reference, overwrite=TRUE) {
    # This warping runs all the time,... while it is fairly fast, it probably can be skipped if all raster layers lign up correctly
    rasters = mapply(
        function(filename, filetype) {
            # Get the nodatavalue if available. Shall set to source's nodatavalue
            gdinfo = rgdal::GDALinfo(filename)
            mdata = attr(gdinfo, 'df')
            dtype = as.character(mdata[['GDType']])
            hasNoDataValues = mdata[['hasNoDataValue']]

            r = bccvl.raster.load(filename)
            # warp, crop and rescale raster file if necessary
            dir = dirname(filename)
            tmpf = file.path(dir, paste0(basename(tempfile()), '.tif')) # TODO: better filename and location?
            te = extent(reference)

            # This is to fix issue with NA value being treated as value 0 if nodatavalue is not set.
            if (hasNoDataValues) {
                # set nodatavalue to the original nodatavalue in the source file
                gdalwarp(filename, tmpf,
                         s_srs=CRSargs(crs(r)), t_srs=CRSargs(crs(reference)),
                         te=c(te@xmin, te@ymin, te@xmax, te@ymax),
                         ts=c(ncol(reference), nrow(reference)),
                         # tr=c(...), ... either this or ts
                         r="near",
                         of="GTiff",
                         dstnodata=mdata[['NoDataValue']],
                         co=c("TILED=YES", "COMPRESS=LZW")
                         )
            }
            else {
                # call gdalwarp without dstnodata
                gdalwarp(filename, tmpf,
                         s_srs=CRSargs(crs(r)), t_srs=CRSargs(crs(reference)),
                         te=c(te@xmin, te@ymin, te@xmax, te@ymax),
                         ts=c(ncol(reference), nrow(reference)),
                         # tr=c(...), ... either this or ts
                         r="near",
                         of="GTiff",
                         co=c("TILED=YES", "COMPRESS=LZW")
                         )
            }

            # put new file back into place
            rasterfilename = tmpf
            if (overwrite) {
                file.rename(tmpf, filename)
                rasterfilename = filename
            }
            # load new file and convert to categorical if required
            r = raster(rasterfilename)
            if (filetype == "categorical") {
                # convert to factor if categorical
                r = as.factor(r)
            }
            return(r)
        },
        raster.filenames, raster.types)
    return(rasters)
}

# raster.filenames : a vector of filenames that will be loaded as rasters
# resamplingflag: a flag to determine which resampling approach to take
bccvl.rasters.to.common.extent.and.resampled.resolution <- function(raster.filenames, raster.types, resamplingflag, selected_layers=NULL, overwrite=TRUE)
{
    # Load rasters and assign CRS if missing
    rasters = lapply(raster.filenames, bccvl.raster.load)

    # determine common raster shape
    reference = bccvl.rasters.common.reference(rasters, resamplingflag, selected_layers)

    # adjust rasters spatially and convert categorical rasters to factors
    rasters = bccvl.rasters.warp(raster.filenames, raster.types, reference, overwrite)

    return(rasters)
}

# return a RasterStack of given vector of input files
# intersecting extent
# lowest or highest resolution depending upon flag
# selected_layers is a list of layers to be considered when determine the resolution of the raster. If none, consider all layers.
bccvl.enviro.stack <- function(filenames, types, layernames, resamplingflag, selected_layers=NULL) {
    # adjust rasters to same projection, resolution and extent
    rasters = bccvl.rasters.to.common.extent.and.resampled.resolution(filenames, types, resamplingflag, selected_layers)
    # stack rasters
    rasterstack = stack(rasters)
    # assign predefined variable names
    names(rasterstack) = unlist(layernames)
    return(rasterstack)
}

# Remove raster object and its associated raster files (i.e. grd and gri) if any
bccvl.remove.rasterObject <- function(rasterObject) {
    raster_filenames = raster_to_filenames(rasterObject, unique = TRUE)
    for (fname in raster_filenames) {
        if (extension(fname)  == '.grd') {
            file.remove(fname, extension(fname, '.gri'))
        }
    }
    rm(rasterObject)
}

bccvl.sp.transform <- function(data, climate.data)
{
    sp <- SpatialPoints(data)
    if (is.na(crs(sp))) {
        crs(sp) <- '+init=epsg:4326'
    }

    # project to the same crs as climate data
    if (!compareCRS(sp, climate.data, verbatim=TRUE)) {
        sp <- spTransform(sp, crs(climate.data))
    }
    return(sp)
}


bccvl.mask <- function(raster, parsedgeojson) {
    # Crop the raster to the extent of the constraint region before masking
    cropped_raster = crop(raster, extent(parsedgeojson))

    # save the constrained raster in work directory instead of raster temporary directory as
    # predict.R clears the raster temp files.
    envraster_filename = paste(bccvl.env$workdir, basename(tempfile(fileext = ".grd")), sep="/")
    masked_raster = mask(cropped_raster, parsedgeojson, filename = envraster_filename)

    # Adjust the levels as some levels may be dropped
    if (is.factor(masked_raster)) {
        masked_raster = as.factor(masked_raster)
    }

    # Remove cropped raster and associated raster files (i.e. grd and gri)
    bccvl.remove.rasterObject(stack(cropped_raster))

    return(masked_raster)
}


# geographically constrained modelling
bccvl.sdm.geoconstrained <- function(rasterstack, occur, absen, rawgeojson, generateCHull) {

    if (is.null(rawgeojson) & !generateCHull) {
        return(list("raster" = rasterstack, "occur" = occur, "absen" = absen))
    }

    # create a dummy geojson for convex-hull polygon if no geojson
    geojson_crs = CRS("+init=epsg:3857")
    if (is.null(rawgeojson)) {
        parsedgeojson <- SpatialPolygons(list(Polygons(list(Polygon(rbind(c(1,1)))), ID=1)), proj4string=crs(rasterstack))
    } else {
        # Parse the geojson from text to SpatialPointsDataFrame
        parsedgeojson <- readOGR(dsn = rawgeojson, layer = "OGRGeoJSON", verbose = FALSE)
        geojson_crs <- crs(parsedgeojson)
    }

    # Assign the same projection to the raster
    if (!compareCRS(rasterstack, parsedgeojson, verbatim=TRUE)) {
        # CRS is different, reproject geojson to rasterstack
        parsedgeojson <- spTransform(parsedgeojson, crs(rasterstack))
    }

    # If there are occurrence points, constraint them
    if (!is.null(occur)) {
        # Constrain the occurrence points
        occurSP <- SpatialPoints(occur)
        # We have to make sure occurSP has the same CRS
        if (is.na(crs(occurSP))) {
            crs(occurSP) <- '+init=epsg:4326'
        }

        if (!compareCRS(occurSP, parsedgeojson, verbatim=TRUE)) {
            occurSP <- spTransform(occurSP, crs(parsedgeojson))
        }

        # actual constraint is the intersection between the occurrence's convex-hull polygon and the constraint.
        # Otherwise, actual constraint is the convex-hull polygon.
        region_offset = 0
        if (generateCHull) {
            # get the offset 
            constraintjson <- rjson::fromJSON(rawgeojson)
            region_offset <- constraintjson$properties$region_offset
            if (is.null(region_offset) || is.na(region_offset) || region_offset == '') {
                region_offset = 0
            }
            else {
                region_offset <- as.double(region_offset)
                region_offset <- ifelse(!is.na(region_offset) && is.numeric(region_offset), region_offset/111.0, 0) # convert from km to degree
            }

            chcoords <- occurSP@coords[chull(occurSP@coords[,1:2]),]
            chullPolygon <- SpatialPolygons(list(Polygons(list(Polygon(chcoords[,1:2], hole=FALSE)), ID=1)), proj4string=crs(parsedgeojson))
            if (!is.null(rawgeojson)) {
                parsedgeojson <- intersect(parsedgeojson, chullPolygon)
            }
            else {
                parsedgeojson <- chullPolygon
            }
        }

        # Add a small buffer of width 1-resolution cell. This is to fix the issue
        # with missing env values along the boundary of the polygon.
        parsedgeojson <- gBuffer(parsedgeojson, byid=TRUE, width=max(region_offset, max(res(rasterstack@layers[[1]]))))

        # Save the convex-hull generated as geojson.
        if (generateCHull) {
            filename = file.path(bccvl.env$outputdir, 'modelling_region.json')
            transformed_geojson = spTransform(parsedgeojson, geojson_crs)
            writeOGR(transformed_geojson, filename, 'OGRGeoJSON', driver='GeoJSON')
            transformed_geojson = NULL

            # Add in the CRS and properties. A bug in writeOGR does not write CRS.
            gjson = rjson::fromJSON(file=filename)
            origjson <- rjson::fromJSON(rawgeojson)
            if (! 'crs' %in% names(gjson)) {
                gjson = append(gjson, list(crs=origjson$crs))   
            }
            if (! 'properties' %in% names(gjson)) {
                gjson = append(gjson, list(properties=origjson$properties))
            }
            write(rjson::toJSON(gjson), filename)

            # clear them
            origjson = NULL
            gjson = NULL
        }    

        fid = names(parsedgeojson)
        cat("\nfid used = ", fid, "\n")
        occurSPconstrained <- occurSP[!is.na(over(occurSP, parsedgeojson)[fid[1]])]
        occurconstrained <- as.data.frame(occurSPconstrained)
        # rest of scripts expects names "lon", "lat" and not "x", "y"
        #names(occurconstrained) <- c("lon", "lat")

        # constraint the true absence points if available
        absenconstrained = NULL
        # Ensure true absence dataset is in same projection system as climate 1st.
        if (!is.null(absen) && nrow(absen) > 0) {
            absenSP <- bccvl.sp.transform(absen, rasterstack)
            absenSPconstrained <- absenSP[!is.na(over(absenSP, parsedgeojson))]

            # project it back to epsg:4326 for saving as a csv file
            absenSPconstrained <- spTransform(absenSPconstrained, CRS('+init=epsg:4326'))
            absenconstrained <- as.data.frame(absenSPconstrained)
            #names(absenconstrained) <- c("lon", "lat")
            #write.csv(absen, file=absenFilename, row.names=FALSE)
        }
    }
    else {
        occurconstrained = NULL
        absenconstrained = NULL
    }

    # Mask the rasterstack (and make sure it is a RasterStack)
    geoconstrained <- stack(lapply(as.list(rasterstack), bccvl.mask, parsedgeojson))

    # Return the masked raster stack and constrained occurrence points
    mylist <- list("raster" = geoconstrained, "occur" = occurconstrained, "absen" = absenconstrained)
    return(mylist)
}

# function to plot projection tiff file (with histogram)
bccvl.plotProjection <- function(inputfile, main) {
    ## define the breaks of the color key
    my.at <- seq(0,1.0,by=0.1)
    ## the labels will be placed vertically centered
    my.labs.at <- seq(0,1.0,by=0.25)
    ## define the labels
    my.lab <- seq(0,1.0,by=0.25)
    ## define colors
    my.col <- colorRampPalette(c("grey90","yellow4","green4"))(100)

    # Read in tiff input file as rasterstack and plot it
    require('rasterVis')
    levelplot(stack(raster(inputfile)),
              at=my.at,
              margin=T,
              col.regions=my.col,
              main=main,
              colorkey=list(labels=list(
                labels=my.lab,
                at=my.labs.at)))
}

# function to generate a filename for the specified file type and extension.
bccvl.get_filepath <- function(file_type, projection_name, species, outputdir=bccvl.env$outputdir, filename_ext=NULL, file_ext='tif') {
    if (is.null(filename_ext)) {
        basename = paste(file_type, projection_name, species, sep="_")
    }
    else {
        basename = paste(file_type, projection_name, species, filename_ext, sep="_")
    }
    return(file.path(outputdir, paste(basename, file_ext, sep=".")))
}

# function to save projection as png image
bccvl.saveProjectionImage <- function(inputfile, projection.name, species, species_algo_str, outputdir=bccvl.env$outputdir, filename_ext=NULL) {
    filename = bccvl.get_filepath("proj", projection.name, species_algo_str, outputdir, filename_ext, "png")
    png(filename)
    title = paste(species, projection.name, "projections", sep=" ")
    plot(raster(inputfile), main=title, xlab='longitude', ylab='latitude')
    # TODO: to use levelplot to produce histogram instead of plot.
    #bccvl.plotProjection(inputfile, title)
    dev.off()
}

# function to compute and save occurrence probability change metrics as geotif file
bccvl.generateOccurrenceProbChangeMetric <- function(prob_rasters, outfilename) {
    changeproj <- overlay(prob_rasters[[1]], prob_rasters[[2]], fun=function(r1, r2) { return(r1-r2) })
    writeRaster(changeproj, outfilename, format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
}

# function to compute and save species range change metric as geotif file
bccvl.generateSpeciesRangeChangeMetric <- function(prob_rasters, threshold, outfilename) {
    # return 1 for Blank, 3 for Expansion, 0 for Contraction and 2 for No Change
    rangeChange <- overlay(as.integer(prob_rasters[[1]] >= threshold),
                           as.integer(prob_rasters[[2]] >= threshold),
                           fun=function(fp, cp) { return((2 * fp) + 1 - cp)})
    writeRaster(rangeChange, outfilename, format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)

    # compute the area for each change category
    grid_area <- raster.from.asc(grid.area(asc.from.raster(rangeChange)))
    total_pixels = ncell(na.omit(as.data.frame(rangeChange)))
    chg_summary = as.data.frame(matrix(ncol=3, nrow=4))
    rownames(chg_summary) <- c('Contraction', 'Blank', 'No Change', 'Expansion')
    colnames(chg_summary) <- c('no_grid_cells', '%_grid_cells', 'area_km2')
    for (i in c(0,1,2,3)) {
        no_pixels = length(rangeChange[rangeChange == i])
        chg_summary[i+1,] <- c(
                no_pixels,
                (no_pixels*100.0)/total_pixels,
                sum(grid_area[rangeChange == i])/1000000.0
            )
    }

    # write it to a file.
    outfilename2 = outfilename
    ext = file_ext(outfilename)
    if (!is.null(ext)) {
        pattern = paste0('\\.', ext, '$')
        outfilename2 <- sub(pattern, '', outfilename)
    }
    outfilename2 = paste(outfilename2, 'csv', sep=".")
    write.csv(chg_summary, file=outfilename2, row.names=TRUE)
}

# function to compute and save Centre of Gravity as csv file.
bccvl.generateCentreOfGravityMetric <- function(projfiles, outfilename) {
    future_proj = raster(projfiles[[1]])
    current_proj = raster(projfiles[[2]])
    future_cog = COGravity(future_proj)
    current_cog = COGravity(current_proj)

    # Do not generate CoG if it has NaN value
    if (is.nan(current_cog['COGy']) || is.nan(current_cog['COGx']) || is.nan(future_cog['COGy']) || is.nan(future_cog['COGx'])) {
        return()
    }

    results = as.data.frame(matrix(ncol=5, nrow=3))
    rownames(results) = c('Centre_of_Range', 'Minimum', 'Maximum')
    colnames(results) = c('current_latitude', 'current_longitude', 'future_latitude', 'future_longitude', 'change_in_m')
    results[1,] = distance(current_cog['COGy'], current_cog['COGx'], future_cog['COGy'], future_cog['COGx'])
    results[2,] = distance(min(coordinates(current_proj)[,2]),
                           min(coordinates(current_proj)[,1]),
                           min(coordinates(future_proj)[,2]),
                           min(coordinates(future_proj)[,1])
                          )
    results[3,] = distance(max(coordinates(current_proj)[,2]),
                           max(coordinates(current_proj)[,1]),
                           max(coordinates(future_proj)[,2]),
                           max(coordinates(future_proj)[,1])
                          )
    write.csv(results, file=outfilename)
}


# function to save projection output raster
bccvl.saveModelProjection <- function(model.obj, projection.name, species, species_algo_str, outputdir=bccvl.env$outputdir, filename_ext=NULL) {
    ## save projections under biomod2 compatible name:
    ##  proj_name_species.tif
    ##  only useful for dismo outputs

    filename = bccvl.get_filepath("proj", projection.name, species_algo_str, outputdir, filename_ext, "tif")
    writeRaster(model.obj, filename, format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)

    # TODO: can we merge this bit with bccvl.saveProjection in eval.R ?
    # Save as image as well
    pngfilename = bccvl.get_filepath("proj", projection.name, species_algo_str, outputdir, filename_ext, "png")
    png(pngfilename)
    title = paste(species, projection.name, "projections", sep=" ")
    plot(model.obj, xlab="latitude", ylab="longtitude", main=title)
    dev.off()

    return (filename)
}

# function to save RData in outputdir
bccvl.save <- function(robj, name, outputdir=bccvl.env$outputdir) {
    filename = file.path(outputdir, name)
    save(robj, file=filename)
}

# function to save CSV Data in outputdir
bccvl.write.csv <- function(robj, name, outputdir=bccvl.env$outputdir, rownames=TRUE) {
    filename = file.path(outputdir, name)
    write.csv(robj, file=filename, row.names=rownames)
}

# function to get model object
bccvl.getModelObject <- function(model.file=bccvl.env$inputmodel) {
    return (get(load(file=model.file)))
}

# convert all .gri/.grd found in folder to gtiff
# TODO: extend to handle other grid file formats, e.g. .asc
bccvl.grdtogtiff <- function(folder, algorithm, filename_ext=NULL, noDataValue=NULL) {
    grdfiles <- list.files(path=folder,
                           pattern="^.*\\.grd")
    for (grdfile in grdfiles) {
        # get grid file name without the extension
        # file_path_sans_ext does not work when file has double '.' before extension i.e. filename..grd
        #grdname <- file_path_sans_ext(grdfile)
        ext = file_ext(grdfile)
        if (!is.null(ext)) {
            pattern = paste0('\\.', ext, '$')
            grdname <- sub(pattern, '', grdfile)
        }

        # read grid raster
        grd <- raster(file.path(folder, grdfile))

        if (is.na(proj4string(grd))) {
            # Projection is missing, initialise it to EPSG:4326
            crs = CRS("+init=epsg:4326")
            proj4string(grd) <- crs
        }

        # write raster as geotiff
        basename = paste(grdname, algorithm, sep="_")
        if (!is.null(filename_ext)) {
            basename = paste(grdname, algorithm, filename_ext, sep="_")
        }
        filename = file.path(folder, paste(basename, 'tif', sep="."))

        # To do: This is a temporary fix for nodatavalue is not recognised by mosaic_raster
        # due to a bug in gdal libarry. It shall be removed when using gdal 2.1.3.
        dtype = dataType(grd)
        if (is.null(noDataValue)) {
            writeRaster(grd, filename, datatype=dataType(grd),
                        format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
        }
        else {
            writeRaster(grd, filename, datatype=dataType(grd), NAflag=noDataValue,
                        format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
        }
        # remove grd files
        file.remove(file.path(folder, paste(grdname, c("grd","gri"), sep=".")))
    }
}

############################################################
#
# define helper functions for projections
#
############################################################

# function to check that the environmental layers used to project the
# model are the same as the ones used to create the model object
#    model.obj     ... model to project
#    climatelayers ... climate data to project onto
bccvl.checkModelLayers <- function(model.obj, climatelayers, climate_filenames) {
    message("Checking environmental layers used for projection")
    # get the names of the environmental layers from the original model
    if (inherits(model.obj, "DistModel")) {
        # dismo package
        model.layers = colnames(model.obj@presence)
    } else if (inherits(model.obj, "gbm")) {
        # brt package
        model.layers = summary(model.obj)$var
    } else if (inherits(model.obj, "BIOMOD.models.out")) {
        # biomod package
        model.layers = model.obj@expl.var.names
    }

    # get the names of the climate scenario's env layers
    pred.layers = names(climatelayers)

    # check if the env layers were in the original model
    if (sum(!(pred.layers %in% model.layers)) > 0 ){
        # To do: Shall remove this sometimes later.
        # The model layer name is to be matched to climate layer name, or its file name.
        # This is to allow for old SDM result to be used.
        if (sum(!(model.layers %in% pred.layers)) > 0){
            filenames = lapply(climate_filenames, function(x) sub("^([^.]*).*", "\\1", basename(x)))
            indexes = match(model.layers, filenames)
            for (i in indexes){
                if (!is.na(i)){
                    pred.layers[i] = model.layers[i]    #Use the corresponding layer name in the model
                }
            }
            names(climatelayers) = pred.layers
        }

        message("Dropping environmental layers not used in the original model creation...")
        # create a new list of env predictors by dropping layers not in the original model
        new.predictors = climatelayers
        for (pl in pred.layers) {
            if (!(pl %in% model.layers)) {
                new.predictors = dropLayer(new.predictors, pl)
            }
        }
        return(new.predictors)
    } else {
        return(climatelayers)
    }
}


family_from_string <- function(s)
{
    # get family from a string (character) in a safe way
    # works for all variants of the R family object (e.g. see ?family)
    # i.e.
    # family_from_string("binomial")
    # family_from_string("binomial(link=logit)")
    # family_from_string("binomial(link=\"logit\")")
    # ...
    # family_from_string("quasi(link = \"identity\", variance = \"constant\")")

    s=gsub(pattern="\"|| ", replacement="", s) # strip quotes and spaces
    f=gsub(pattern="\\(.*\\)", replacement="", s) # the name of the function

    allowable= c("binomial",
                "gaussian",
                "Gamma",
                "inverse.gaussian",
                "poisson",
                "quasi",
                "quasibinomial",
                "quasipoisson")

    if (! f %in% allowable )
    {
        stop(sprintf("unsupported function %s", f))
    }

    fargs=gsub(pattern=".*\\(||\\)",
               replacement="",
               sub(pattern=f,
                    replacement="",
                    s)) #get the args inside the parentheses
    args=list()

    if (fargs != "")
    {
        l=strsplit(fargs, ",")[[1]]
        for( i in 1:length(l) )
        {
            ll=strsplit(l[i],"=")[[1]]
            if (length(ll) == 2)
            {
                args[ll[1]] = ll[2]
            }
            else
            {
                stop(sprintf("unhandled result when splitting %s", l[i]))
            }
        }
    }
    return (do.call(what=f, args=args))
}

#' Grid Information from Geographic (lat lon) Projections
#'
#' Since spatial grids in geographic projections do not have equal area or
#' perimeters, \code{grid.info} extracts perimeter & area related information
#' for latitudinal bands with differing longitudinal widths. \cr\cr Outputs
#' lengths are in m using Vincenty's equation (\code{distance})and areas in m2.
#' Surface areas are calculated summing surface areas of spherical polygons as
#' estimated using l'Huiller's formula.
#'
#'
#' @param lats is a vector of latitudes representing the midpoint of grid cells
#' @param cellsize is a single value (assuming square cells) or a two value
#' vector (rectangular cells) representing the height (latitude) and width
#' (longitude) of the cells
#' @param r is a single value representing the radius of the globe in m.
#' Default is for the WGS84 elipsoid
#' @return a data.frame listing: \item{lat}{the latitude representing the
#' midpoint of the cell} \item{top}{length of the top of the cell (m)}
#' \item{bottom}{length of the bottom of the cell (m)} \item{side}{length of
#' the side of the cell (m)} \item{diagnal}{length of the diagnals of the cell
#' (m)} \item{area}{area of the cell (m2)}
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @references information on l'Huiller's formula
#' \url{http://williams.best.vwh.net/avform.htm for more info)} code for
#' estimating area of polygon on sphere was modified from
#' \url{http://forum.worldwindcentral.com/showthread.php?t=20724}
#' @examples
#'
#' #show output for latitudes from -87.5 to 87.5 at 5 degree intervals
#' grid.info(lats=seq(-87.5,87.5,5), 5)
#'
#' @export
# This is a fix for SDM tool grid.info function due to floating point operation.
grid.info <- function(lats,cellsize,r=6378137) {
    r2 = r^2 #radius of earth
    ###need checks to ensure lats will not go beyond 90 & -90
    if (length(cellsize)==1) cellsize=rep(cellsize,2) #ensure cellsize is defined for both lat & lon
    out = data.frame(lat=lats) #setup the output dataframe
    toplats = lats+(0.5*cellsize[1]); bottomlats = lats-(0.5*cellsize[1]) #define the top and bottom lats
    check = range(c(toplats,bottomlats),na.rm=TRUE); if (-90.0001>check[1] | 90.0001<check[2]) stop('latitudes must be between -90 & 90 inclusively')
    out$top = distance(toplats,rep(0,length(lats)),toplats,rep(cellsize[2],length(lats)))$distance
    out$bottom = distance(bottomlats,rep(0,length(lats)),bottomlats,rep(cellsize[2],length(lats)))$distance
    out$side = distance(toplats,rep(0,length(lats)),bottomlats,rep(0,length(lats)))$distance
    out$diagnal = distance(toplats,rep(0,length(lats)),bottomlats,rep(cellsize[2],length(lats)))$distance
    #calculate area of a spherical triangle using spherical excess associated by knowing distances
    #tan(E/4) = sqrt(tan(s/2)*tan((s-a)/2)*tan((s-b)/2)*tan((s-c)/2))
    #where a, b, c = sides of spherical triangle
    #s = (a + b + c)/2
    #from CRC Standard Mathematical Tables
    #calculate excess based on  l'Huiller's formula (http://williams.best.vwh.net/avform.htm for more info)
    #code modified from (http://forum.worldwindcentral.com/showthread.php?t=20724)
    excess = function(lam1,lam2,beta1,beta2){ #calculate excess... inputs are in radians
        haversine = function(y) { (1-cos(y))/2 }
        cosB1 = cos(beta1); cosB2 = cos(beta2)
        hav1 = haversine(beta2-beta1) + cosB1*cosB2*haversine(lam2-lam1)
        aa = 2 * asin(sqrt(hav1)); bb = 0.5*pi - beta2; cc = 0.5*pi - beta1
        ss = 0.5*(aa+bb+cc)
        tt = tan(ss/2)*tan((ss-aa)/2)*tan((ss-bb)/2)*tan((ss-cc)/2)
        return(abs(4*atan(sqrt(abs(tt)))))
    }
    if (any(bottomlats==-90)) { pos = which(bottomlats==-90); bottomlats[pos] = -bottomlats[pos]; toplats[pos] = -toplats[pos]} #ensure no -90 bottom lats
    out$area = excess(lam1=0,lam2=cellsize[2]*pi/180,toplats*pi/180,toplats*pi/180)
    out$area = abs(out$area-excess(lam1=0,lam2=cellsize[2]*pi/180,bottomlats*pi/180,bottomlats*pi/180))*r2
    return(out)
}

#' Vincenty Direct Calculation of Distance and Direction
#'
#' \code{distance} estimates the distance given a starting & ending latitude
#' and longitude. \cr \cr For general information on Vincenty's formula, see
#' e.g., \url{http://en.wikipedia.org/wiki/Vincenty's_formulae}. It states: \cr
#' \emph{Vincenty's formulae are two related iterative methods used in geodesy
#' to calculate the distance between two points on the surface of an spheroid,
#' developed by Thaddeus Vincenty in 1975. They are based on the assumption
#' that the figure of the Earth is an oblate spheroid, and hence are more
#' accurate than methods such as great-circle distance which assume a spherical
#' Earth.} \cr \cr \bold{Note:} this method assumes a locations are lat & lon
#' given in WGS 84.\cr\cr Direction, if requested, is the the initial bearing
#' (sometimes referred to as forward azimuth) for which one would follow as a
#' straight line along a great-circle arc from start to finish.\cr \cr
#' \bold{Note:} this will fail if there are NA's in the data.
#'
#'
#' @param lat1 a single value or vector of values representing latitude in
#' decimal degrees from -90 to 90 degrees. Alternatively, a data.frame or
#' matrix can be used here with each column representing lat1, lon1, lat2, lon2
#' (in that order).
#' @param lon1 a single value or vector of values representing longitude in
#' decimal degrees from -180 to 180 degrees. If NULL, lat1 is assumed to be a
#' matrix or data.frame.
#' @param lat2 a single value or vector of values representing latitude in
#' decimal degrees from -90 to 90 degrees. If NULL, lat1 is assumed to be a
#' matrix or data.frame.
#' @param lon2 a single value or vector of values representing longitude in
#' decimal degrees from -180 to 180 degrees. If NULL, lat1 is assumed to be a
#' matrix or data.frame.
#' @param bearing boolean value as to calculate the direction as well as the
#' distance.
#' @return Returns a data.frame with: \item{lon1}{the original longitude}
#' \item{lat1}{the original latitude} \item{lon2}{the destination longitude}
#' \item{lat2}{the destination latitude} \item{distance}{the distance used}
#' \item{bearing}{if requested, the bearing between the two points}
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @seealso \code{\link{destination}}
#' @references Vincenty, T. 1975. Direct and Inverse Solutions of Geodesics on
#' the Ellipsoid with application of Nested Equations. Survey Review, vol XXII
#' no 176. \url{http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf}
#' @source The source code for the distance algorithm here was modified from
#' \url{http://www.movable-type.co.uk/scripts/latlong-vincenty.html}.\cr \cr
#' Distances were validated against Geoscience Australia calculations
#' (\url{http://www.ga.gov.au/geodesy/datums/vincenty_inverse.jsp}).\cr \cr
#' Bearings were from multiple sources including
#' \url{http://williams.best.vwh.net/avform.htm#Crs}.
#' @examples
#'
#'
#' #get the distance of 1 degree longitude at each 5 degrees latitude from -90 to 90
#' distance(lat1=seq(-90,90,5),lon1=rep(0,37),lat2=seq(-90,90,5),lon2=rep(1,37),bearing=TRUE)
#'
#'
#' @export
#' @useDynLib SDMTools Dist
# This is a fix for SDM tool distance function due to floating point operation.
distance = function(lat1, lon1=NULL, lat2=NULL, lon2=NULL, bearing=FALSE) {
    if (is.data.frame(lat1) | is.matrix(lat1)) { #if input is matrix or data.frame... break it out to individual vectors
        lat1 = as.matrix(lat1); if (ncol(lat1)!=4) stop('incorrect lat/lon inputs... must be matrix with 4 columns or 4 vectors')
        lon2=lat1[,4]; lat2=lat1[,3]; lon1=lat1[,2]; lat1=lat1[,1] #break out individual columns
    } else if (!is.null(lat2) & !is.null(lon1) & !is.null(lon2)) {
        if (!all(c(length(lat2),length(lon1),length(lon2))==length(lat1))) stop('inputs must all be of same length')
    } else { stop('inappropriate inputs... see helpfile') }
    if (any(c(lon1,lon2) < -180.0001) | any(c(lon1,lon2) > 180.0001)) stop('lon must be decimal degrees between -180 & 180')
    if (any(c(lat1,lat2) < -90.0001) | any(c(lat1,lat2) > 90.0001)) stop('lat must be decimal degrees between -90 & 90')
    #cycle through and output the new data
    out = data.frame(lat1=lat1,lon1=lon1,lat2=lat2,lon2=lon2)
    out$distance = round(.Call('Dist',out$lat1,out$lon1,out$lat2,out$lon2,PACKAGE='SDMTools'),2) #round to the nearest mm
    if (bearing) { #if requested, calculate bearing
        lat1=lat1*pi/180;lat2=lat2*pi/180;lon1=lon1*pi/180;lon2=lon2*pi/180 #convert to radians
        brng = atan2(sin(lon2-lon1)*cos(lat2),cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(lon1-lon2)) #estimate bearing
        out$bearing = ((brng*180/pi)+360)%%360 #convert to bearing in degrees
    }
    #return the output
    return(out)
}


########################################################
# Patch Biomod2 sample.factor.levels function
#  taken from: biomod2-3.3-7 : https://github.com/cran/biomod2/blob/master/R/sample.factor.levels.R#L70
########################################################

sample.factor.levels <- function(x, mask.out = NULL, mask.in = NULL){
  ## make some checking of given parameters
  ## TODO(damien)
  if(inherits(x, 'Raster')){
    fact.level.cells <- .sample.factor.levels.raster(x, mask.out = mask.out, mask.in = mask.in)
    return(fact.level.cells)
  } else if(inherits(x, 'data.frame')){
    fact.level.cells <- .sample.factor.levels.data.frame(x, mask.out = mask.out, mask.in = mask.in)
    return(fact.level.cells)
  } else {
    warning(paste0("\nunsupported input data.",
                   "\nx should be a Raster* object or a data.frame.",
                   "\n NULL returned"))
    return(NULL)
  }
}

.sample.factor.levels.raster = function (x, mask.out = NULL, mask.in = NULL) 
{
  fact.var <- which(is.factor(x))
  if (any(fact.var)) {
    fact.level.cells <- sapply(fact.var, function(f) {
      selected.cells <- NULL
      fact.level.original <- unlist(raster::levels(subset(x, 
                                                          f)))
      fact.level <- fact.level.original
      cat("\n> fact.level for", names(x)[f], ":\t", paste(fact.level, 
                                                          names(fact.level), sep = ":", collapse = "\t"))
      if (!is.null(mask.out)) {
        fact.levels.sampled <- unlist(levels(as.factor(mask(subset(x, 
                                                                   f), mask.out))))
        attr(fact.levels.sampled, "names") <- attr(fact.level.original, 
                                                   "names")[fact.levels.sampled]
        cat("\n - according to mask.out levels", fact.levels.sampled, 
            "have already been sampled")
        fact.level <- fact.level[!is.element(fact.level,
                                             fact.levels.sampled)]
      }
      if (length(fact.level)) {
        if (!is.null(mask.in)) {
          for (mask.in.id in 1:length(mask.in)) {
            if (length(fact.level)) {
              x.f.masked <- as.factor(mask(subset(x, 
                                                  f), mask.in[[mask.in.id]]))
              x.f.levels <- unlist(levels(x.f.masked))
              attr(x.f.levels, "names") <- attr(fact.level.original, 
                                                "names")[x.f.levels]
              fact.levels.in.m.in <- fact.level[is.element(fact.level, 
                                                           x.f.levels)]
              if (length(fact.levels.in.m.in)) {
                cat("\n - levels", fact.levels.in.m.in, 
                    "will be sampled in mask.out", mask.in.id)
                selected.cells <- c(selected.cells, sapply(fact.levels.in.m.in, 
                                                           function(fl) {
                                                             sample(which(x.f.masked[] == fl), 
                                                                    1)
                                                           }))
                fact.level <- fact.level[!is.element(fact.level, 
                                                     fact.levels.in.m.in)]
              }
            }
          }
        }
        if (length(fact.level)) {
          cat("\n - levels", fact.level, "will be sampled in the original raster")
          selected.cells <- c(selected.cells, sapply(fact.level, 
                                                     function(fl) {
                                                       sample(which(subset(x, f)[] == fl), 1)
                                                     }))
        }
      }
      return(selected.cells)
    })
    fact.level.cells <- as.numeric(fact.level.cells[-which(sapply(fact.level.cells, is.null))])
    return(fact.level.cells)
  }
  else {
    return(NULL)
  }
}

.sample.factor.levels.data.frame <- function(x, mask.out = NULL, mask.in = NULL){
  ## identificate the factorial variables
  fact.var <- which(sapply(x, is.factor))
  ## check if some factorial variables are in the input data
  if(any(fact.var)){ ## some factorial variables present
    fact.level.cells <- as.numeric(unlist(sapply(fact.var, function(f){
      ## initialize the list of cells that are selected
      selected.cells <- NULL
      ## get the levels of the factor on the full dataset
      fact.level.original <- levels(x[, f])
      fact.level <- fact.level.original
      cat("\n> fact.level for",  colnames(x)[f], ":\t", paste(1:length(fact.level), fact.level, sep = ":", collapse = "\t"))
      if(!is.null(mask.out)){ ## mask containing points that have already been sampled
        ## check the levels of the fector that have been already sampled
        fact.levels.sampled <- unique(na.omit(as.character(x[mask.out[, 1], f])))
        ## remove already sampled points from candidates
        x[mask.out[, 1], ] <- NA
        #         ## update levels names (lost during mask conversion)
        #         attr(fact.levels.sampled, "names") <- attr(fact.level.original, "names")[fact.levels.sampled]
        cat("\n - according to mask.out levels", fact.levels.sampled, "have already been sampled")
        ## update the list of factor levels to sample
        fact.level <- setdiff(fact.level, fact.levels.sampled)
      }
      if(length(fact.level)){
        ## try first to sample factors in the given masks
        if(!is.null(mask.in)){ ## list of mask we want to sample in (order matter!)
          for(mask.in.id in 1:ncol(mask.in)){
            ## check that some levels remains to be sampled
            if(length(fact.level)){
              ## update the masked version of the factorial raster
              x.f.masked <- as.character(x[, f])
              x.f.masked[!mask.in[, mask.in.id]] <- NA
              x.f.levels <- unique(na.omit(x.f.masked))
              ## get the list of levels that coulb be sampled in this mask
              fact.levels.in.m.in <- intersect(fact.level, x.f.levels)
              if(length(fact.levels.in.m.in)){
                cat("\n - levels", fact.levels.in.m.in, "will be sampled in mask.out", mask.in.id)
                selected.cells <- c(selected.cells, sapply(fact.levels.in.m.in, function(fl){
                  candidate.cells <- na.omit(which(x.f.masked[] == fl))
                  selected.cell <- NULL
                  if(length(candidate.cells) == 1){ ## single candiate cell
                    selected.cell <- candidate.cells
                  } else { ## multi candidate cells
                    selected.cell <- sample(candidate.cells, 1)
                  }
                  return(selected.cell)
                }))
                ## update the list of factor levels to sample
                fact.level <- setdiff(fact.level, fact.levels.in.m.in)
              }
            } 
          } ## end loop over mask.in
        }
        ## @note if some levels remains unsampled then we will take a random value of
        ## them in the full dataset => !! this should be tricky if mask.in arg is given
        ## because the value will be picked out of mask.in but is necessary to 
        ## ensure that models will run smoothly
        if(length(fact.level)){
          cat("\n - levels", fact.level, "will be sampled in the original data.frame")
          selected.cells <- c(selected.cells, sapply(fact.level, function(fl){
            candidate.cells <- na.omit(which(x[, f] == fl))
            selected.cell <- NULL
            if(length(candidate.cells) == 1){ ## single candiate cell
              selected.cell <- candidate.cells
            } else { ## multi candidate cells
              selected.cell <- sample(candidate.cells, 1)
            }
            return(selected.cell)
          }))
        }
      }
      return(selected.cells)
    })))
    return(unique(fact.level.cells))
  } else { ## no factorial variable
    return(NULL)
  }
}
assignInNamespace("sample.factor.levels",sample.factor.levels, ns="biomod2")
assignInNamespace("grid.info",grid.info, ns="SDMTools")
assignInNamespace("distance",distance, ns="SDMTools")

### BCCVL SDM model evaluation script ###

### Includes:
# 1: Functions to save output
# 2: Function to evaluate SDM model (written by A/Prof. Sama Low-Choy, Griffith University)
# 3: Functions to create model outputs
# 4: Run the evaluation and save outputs

#########################################################################
### 1: Functions to save SDM outputs
#########################################################################

bccvl.saveModelEvaluation <- function(out.evaluation, out.stats, out.lossfunction, species_algo_str){
  bccvl.write.csv(data.frame(out.evaluation), 
                  name = paste0(paste("Evaluation-data", species_algo_str, sep="_"), ".csv"))
  bccvl.write.csv(data.frame(out.stats), 
                  name = paste0(paste("Evaluation-statistics", species_algo_str, sep="_"), ".csv"))
  bccvl.write.csv(data.frame(out.lossfunction), 
                  name = paste0(paste("Loss-function-intervals-table", species_algo_str, sep="_"), ".csv"))
  }

bccvl.saveProjection <- function(proj.model, species_algo_str, filename_ext=NULL) {
  if (!is.null(filename_ext)) {
    basename = paste("proj", 'current', species_algo_str, filename_ext, sep="_")
  }
  else {
    basename = paste("proj", 'current', species_algo_str, sep="_")
  }
  png(file=file.path(bccvl.env$outputdir, paste(basename, 'png', sep=".")))
  plot(proj.model, on_0_1000=FALSE)
  dev.off()
}

#########################################################################
### 2: Function to evaluate the SDM model
#########################################################################

absmean <- function(x) abs(mean(x, na.rm=T))
absdiff <- function(x) abs(diff(x, na.rm=T))
pick1min <- function(x) { # w<-c(5,6,7,11,13)
	w <- which(x == min(x, na.rm=T)); 
	len <- length(w); 
	if (len == 1) return(w) else {if (len %% 2 ==1) {return(median(w))} else {return(median(w[-1]))}}
}

performance.2D <- function(obs, pred, species_algo_str, make.plot="bccvl", kill.plot=T) {
  library(gridExtra)
  library(pROC)
  
  # AIM: Calculate 2D measures of predictive performance for any
  # model that predicts a probability of presence (or success), 
  # to be compared to binary observations of presence/absence (or 
  # success/failure).
  #
  # AUTHOR: S.Low-Choy, Jan 2016, Griffith University
  # ACKNOWLEDGEMENTS: Shawn Laffan for useful discussions, 
  # Chantal Huijbers, Sarah Richmond and Linda for testcase
  #
  # INPUTS
  # obs = vector of observations of presence/absence (binary)
  # pred = vector of predicted probabilities of presence
  # 
  # OUTPUTS
  # Predicted presence/absence density plot and histogram
  # Sensitivity/Specificity plot
  # ROC plot
  # Plot with four different error rates across varying treshold probability values
  # 4 different loss functions (maximizing TPR+TNR, balancing all errors, or equalising error rates)
  # Table with best probability threshold value corresponding to minimum loss
  # Range of threshold probability values for which loss falls within 5% of minimum
  # 
  
  # TESTING 23 July 2018 - SLC - problem with threshold = 0 when predictions only 0 or 1
  # species_algo_str <- "SRE"; make.plot <- "bccvl"; kill.plot <- F; obs<-gobs; pred<-gpred
  
  
  # TESTING 23 July 2018 - SLC - when more than one "best" value as the minimum is achieved twice!
  # Addition Number 2
  
  
  #############################################################
  #
  # ERROR CHECKING
  #
  # Check that observations and predictions match length
  nobs <- length(obs)
  if (nobs != length(pred)) stop("Ensure that vectors for observations and predictions are of equal length!")
  
  # Check observations are binary, then recode as a factor if necessary
  tbl.obs <- table(obs)
  diversity.obs <- length(tbl.obs)
  if (diversity.obs==1) stop("All observations are the same!")
  if (diversity.obs!=2) stop("Ensure that the vector of observations only has two possible values.")
  
  # Check coding of presence and absence
  if (is.factor(obs)) {
    truth <- as.character(obs)
    temp <- as.numeric(truth)
    if (!is.na(temp)) 
      truth <- temp 
  } else {
    truth <- obs
  }
  if (is.numeric(truth)) {
    if (any(truth==0)) {
      # presume 0 is coded to be absence
      # if presences are not 1s then recode them to be 1s
      if (any(truth[truth!=0]!=1)) {
        truth[truth!=0] <- 1
      } 
    } else if (any(truth ==1)) {
      # if there are no zeros, then presume 1-2 coding
      if (any(truth ==2)) {
        truth <- truth-1
      }
    } else {
      stop("Can't figure out coding of absence-presence as it is not 0-1 or 1-2. Suggest you recode obs as a factor, with levels 0 and 1")
    } # end if any truth == 0
  } else if (is.character(truth)) {
    # look for "p" for presence, "d" for detected or "o" for occupied
    the.pres <- grep("^[pP]", truth)
    if (any(the.pres)) {
      letter.pres <- "p"
    } else {
      the.pres <- grep("^[Dd]", truth)
      if (any(the.pres)) {
        letter.pres <- "d"
      } else {
        the.pres <- grep("^[oO]", truth)
        if (any(the.pres)) {
          letter.pres <- "o"
        } else {
          stop("Can't figure out coding of presences as they do not start with the letter p, d or o. Suggest you recode presences using one of these options.")
        }
      }
    } # end if any the.pres
    truth[the.pres] <- 1
    
    # recode all non "p" or "d" to be absence (could be "a" for absence, or "n" for not present/seen/detected/occupied, or "u" for unseen etc, or "b" for background)
    tbl.obs <- table(truth)
    letter.abs <- names(tbl.obs)[-grep("1", names(tbl.obs))]
    truth[grep(paste("^",letter.abs,sep=""), truth)] <- 0
    
    if (any(the.pres)) {
      truth <- as.numeric(truth)
    }
  } else {
    stop("Ensure that the data type of observations is numeric, character or factor.")  
  } # end if is.numeric(truth) or is.character(truth)
  
  # Check predictions are probabilities between 0 and 1
  if (any(pred < 0) | any(pred > 1)) stop("Predictions should be probabilities between zero and one. (Check that predictions are not on the log odds scale.)")
  
  #
  # MEASURES
  #
  
  # CREATE ERROR MATRICES
  list.tpv <- sort(unique(pred)) 
  
  tbl.pred <- table(pred)
  diversity.pred <- length(tbl.pred)
  if (diversity.pred==1) stop("All predictions are the same!")
  if (any(list.tpv==0)) { if (diversity.pred==2) { list.tpv <- 1; } else { list.tpv <- list.tpv[list.tpv!=0]; } }
  
  # tpv = threshold probability value: threshold that is used to transform the continuous probability of presence predicted by the model into a binary prediction:
  # probabilities < tpv = absences, probabilities > tpv = presences
  
  tp <- fp <- tn <- fn <- rep(NA, length(list.tpv))
  tpr <- fpr <- tnr <- fnr <- rep(NA, length(list.tpv))
  ppv <- npv <- fdr <- fors <- rep(NA, length(list.tpv))
  acc <- mcr <- tss <- bs <- rep(NA, length(list.tpv))
  tp.rand <- ets <- or <- csi <- rep(NA, length(list.tpv))
  Po <- Pe <- kappa <- roc <- auc <- rep(NA, length(list.tpv))
  
  # CALCULATE 1D MEASURES OF PREDICTIVE PERFORMANCE
  
  ## Elements of contigency table:
  
  # tp = True Positives (observed and predicted presences)
  # fp = False Positives (observed absences predicted as presences)
  # tn = True Negatives (observed and predicted absences)
  # fn = False Negatives (observed presences predicted as absences)
  
  for (ell in seq(along=list.tpv)) {  # ell <- 1
    
    th <- list.tpv[ell]
    
    tp[ell] <- length(which(pred>=th & truth==1))
    fp[ell] <- length(which(pred>=th & truth ==0))
    tn[ell] <- length(which(pred<th & truth ==0))
    fn[ell] <- length(which(pred<th & truth ==1))
    
  ## Evaluation statistics:
    
    # tpr = True Positive Rate (= Sensitivity) = proportion of observed presences that are correctly predicted.
    tpr[ell] <- tp[ell]/(tp[ell]+fn[ell])
    
    # fpr = False Positive Rate = proportion of observed absences that are incorrectly predicted.
    fpr[ell] <- fp[ell]/(fp[ell]+tn[ell])
    
    # tnr = True Negative Rate (= Specificity) = proportion of observed absences that are correctly predicted.
    tnr[ell] <- tn[ell]/(tn[ell]+fp[ell])
    
    # fnr = False Negative Rate = proportion of observed presences that are incorrectly predicted.
    fnr[ell] <- fn[ell]/(fn[ell]+tp[ell])
    
    # Note: tpr + fnr = 1, tnr + fpr = 1
    
    # ppv = Positive Predictive Value = for all predicted presences, how many were true observed presences?
    ppv[ell] <- tp[ell]/(tp[ell]+fp[ell])
    
    # npv = Negative Predictive Value = for all predicted absences, how many were true absences?
    npv[ell] <- tn[ell]/(tn[ell]+fn[ell])
    
    # fdr = False Discovery Rate = for all predicted presences, how many were false presences (observed absences)?
    fdr[ell] <- fp[ell]/(tp[ell]+fp[ell])
    
    # fors = False Omission Rate = for all predicted absences, how many were false absences (observed presences)?
    fors[ell] <- fn[ell]/(tn[ell]+fn[ell])
    
    # Note: ppv + fdr = 1, npv + fors = 1
    
    # acc = Accuracy = proportion of correctly predicted cases.
    acc[ell] <- (tp[ell]+tn[ell])/(tp[ell]+fp[ell]+tn[ell]+fn[ell])
    
    # mcr = Misclassification Rate = proportion of incorrectly predicted cases.
    mcr[ell] <- (fp[ell]+fn[ell])/(tp[ell]+fp[ell]+tn[ell]+fn[ell])
    
    # Note: acc + mcr = 1
    
    # tss = True Skill Statistic - describes how well the model separates presences from absences.
    tss[ell] <- (tpr[ell]-fpr[ell]) 
    
    # bs = Bias Score = frequency of predicted presences compared to the frequency of observed presences.
    bs[ell] <- (tp[ell]+fp[ell])/(tp[ell]+fn[ell])
    
    # csi = Critical Success Index = proportion of observed and predicted presences that are correct.
    csi[ell] <- tp[ell]/(tp[ell]+fp[ell]+fn[ell])
    
    # ets = Equitable Threat Score = proportion of observed and predicted presences that are correct, adjusted for true positives with random chance.
    tp.rand[ell] <- ((tp[ell]+fn[ell])*(tp[ell]+fp[ell]))/(tp[ell]+fp[ell]+tn[ell]+fn[ell])
    ets[ell] <- (tp[ell]-tp.rand[ell])/(tp[ell]+fp[ell]+fn[ell]-tp.rand[ell])  
    
    # or = Odds Ratio = ratio of a correct prediction to an incorrect prediction. 
    or[ell] <- (tp[ell]*tn[ell])/(fp[ell]*fn[ell])
    
    # kappa = Cohen's Kappa = Accuracy of the prediction relative to that of random chance.
    Po[ell] <- acc[ell] # observed accuracy
    Pe[ell] <- ((((tp[ell]+fp[ell])/(tp[ell]+fp[ell]+tn[ell]+fn[ell]))*((tp[ell]+fn[ell])/(tp[ell]+fp[ell]+tn[ell]+fn[ell]))) + (((fn[ell]+tn[ell])/(tp[ell]+fp[ell]+tn[ell]+fn[ell]))*((fp[ell]+tn[ell])/(tp[ell]+fp[ell]+tn[ell]+fn[ell])))) # expected accuracy by random chance
    kappa[ell] <- ((Po[ell]-Pe[ell])/(1-Pe[ell]))
  }

  # auc = Area Under the (ROC) Curve
  roc <- roc(truth, pred)
  auc <- auc(roc)
    
  # Compile the information into dataframes
  temp <- data.frame(list(tpv=list.tpv, tpr=tpr, fpr=fpr,  tnr=tnr, fnr=fnr, ppv=ppv, fdr=fdr, npv=npv, fors=fors))
  auc.d <- data.frame(auc)
  auc.d <- round(auc.d, digits = 2)
  
  # CALCULATE 2D MEASURES OF PREDICTIVE PERFORMANCE
  
  # Calculate losses as cost functions across errors
  list.errors <- c("fpr","fnr","fors","fdr")
  
  # Loss functions:
  # L.diag = minimizing the sum of the diagnostic errors (FPR, FNR)
  #        = maximizing the sum of the True Positive Rate (Sensitivity) and the True Negative Rate (Specificity)
  # L.pred = minimizing the sum of the predictive errors (FDR, FORS)
  #        = maximizing the sum of the Positive Predictive Value and the Negative Predictive Value
  # L.all = balancing all error rates: FPR, FNR, FDR, FORS
  # L.eq.diag = equalising diagnostic errors: FPR = FNR = cross-over of TPR and TNR
  
  temp$L.diag <- apply(temp[, c("fpr","fnr")],1, absmean)
  temp$L.pred <- apply(temp[, c("fors","fdr")],1, absmean)
  temp$L.all <- apply(temp[, list.errors],1, absmean)
  temp$L.eq.diag <- apply(temp[, c("tpr", "tnr")], 1, absdiff)

  # Addition Number 2, 23 July 2018
  # Check if there is more than one minimum, then pick the middle 
  best <- list(
    diag = temp$tpv[pick1min(temp$L.diag) ],
    pred = temp$tpv[pick1min(temp$L.pred == min(temp$L.pred, na.rm=T))],
    all = temp$tpv[pick1min(temp$L.all == min(temp$L.all, na.rm=T))],   
    eq.diag = temp$tpv[pick1min(temp$L.eq.diag==min(temp$L.eq.diag, na.rm=T))]
  )
  
  # End Addition Number 2, 23 July 2018
  
  # Calculate the range of threshold probability values for which each of the losses fall within 5% of the best value
  rangeperf <- matrix(NA, nrow=4, ncol=2, dimnames=list(names(best), c("lower","upper")))
  for (v in names(best)) { # v<-"eq.diag"
    the.v <- paste("L.", v, sep="")
    min.v <- min(temp[,the.v], na.rm=T)
    d.minv <- (abs(temp[,the.v] - min.v) / min.v) 
    the.range <- temp$tpv[ d.minv < 0.05 ]
    rangeperf[dimnames(rangeperf)[[1]]==v, 1:2] <- c(min(the.range, na.rm=T), max(the.range, na.rm=T))
  }
  
  rangeperf <- as.data.frame(rangeperf)
  rangeperf$type.of.loss <- names(best)
  rangeperf$best <- unlist(best)
  
  loss.table <- subset(rangeperf, select = c("lower", "upper", "best"))
  row.names(loss.table) = c("Maximize TPR+TNR", "Maximize PPV+NPV", "Balance all errors", "TPR = TNR")
  loss.table <- as.data.frame(loss.table)
  
  # Rescale
  temp$L.eq.diag <- temp$L.eq.diag/max(temp$L.eq.diag, na.rm=T) 
  
  #########################################################################
  ### 3: Functions to create SDM outputs
  #########################################################################
  
  if (make.plot!="") {
    library(reshape2)
    # reshape the data so that it is in long rather than wide format (= each row represents one item, labels are specified by 'measure' column; used by ggplot2)
    errs <- melt(temp, id.var="tpv", measure.var=c("tpr", "tnr", "fpr", "fnr", "fdr", "fors", "L.diag", "L.pred", "L.all", "L.eq.diag"))
    names(errs)[2] <- c("measure")
    
    # Create Presence/absence density plot across threshold probability values
    temp2 <- data.frame(list(pred=pred, obs=obs))
    png(file=file.path(bccvl.env$outputdir, sprintf("%s-presence-absence-plot_%s.png", make.plot, species_algo_str)), width=480, height=480)
    g1 <- ggplot(temp2, aes(x=pred, fill=factor(obs))) + 
      geom_density(stat="density", alpha=0.5) + 
      labs(title="Presence/absence density plot \nacross predicted probability of presence", x="\nPredicted probability of presence", y="Density\n") + 
      scale_fill_manual(values=c("#EE3B3B", "#6495ED"), labels=c(" Absences      ", " Presences")) +
      theme(axis.text = element_text(family="Arial", size=rel(1.5)), axis.title = element_text(family="Arial", size=rel(1.5)), plot.title = element_text(family="Arial", size=rel(2)), legend.text = element_text(family="Arial", size=rel(1.5)), legend.position="top", legend.key=element_blank(), legend.key.size=unit(1.5, "lines")) + 
      guides(fill=guide_legend(nrow=1, title=NULL))
    print(g1)
    dev.off()
    
    # Create Presence/absence histogram across threshold probability values
    png(file=file.path(bccvl.env$outputdir, sprintf("%s-presence-absence-hist_%s.png", make.plot, species_algo_str)), width=480, height=480)
    g2 <- ggplot(temp2, aes(x=pred, fill=factor(obs)))  + 
      geom_histogram(position="dodge", alpha = 0.5) +
      labs(title="Presence/absence histogram \nacross predicted probability of presence", x="\nPredicted probability of presence", y="Count\n") +
      scale_fill_manual(values=c("#EE3B3B", "#6495ED"), labels=c(" Absences    ", " Presences")) +
      theme(axis.text = element_text(family="Arial", size=rel(1.5)), axis.title = element_text(family="Arial", size=rel(1.5)), plot.title = element_text(family="Arial", size=rel(2)), legend.text = element_text(family="Arial", size=rel(1.5)), legend.position="top", legend.key=element_blank(), legend.key.size=unit(1.5, "lines")) + 
      guides(fill=guide_legend(nrow=1, title=NULL))
    print(g2)
    dev.off()
    
    # Create TPR-TNR plot
    png(file=file.path(bccvl.env$outputdir, sprintf("%s-TPR-TNR_%s.png", make.plot, species_algo_str)), width=480, height=480)
    g3 <- ggplot(errs[errs$measure %in% c("tpr", "tnr"), ], 
                 aes(x=tpv, y=value, colour=measure)) + 
      geom_line(size=1.2) + 
      ylim(0,1) +
      labs(title="Sensitivity-Specificity plot\n", x="\nThreshold probability value", y="TPR/TNR value\n") +
      scale_colour_manual(values=c("#3CAB34", "#049CE3"), labels=c("True Positive Rate (=Sensitivity)", "True Negative Rate (=Specificity)")) + 
      theme(axis.text = element_text(family="Arial", size=rel(1.5)), axis.title = element_text(family="Arial", size=rel(1.5)), plot.title = element_text(family="Arial", size=rel(2)), legend.text = element_text(family="Arial", size=rel(1.5)), legend.position="top", legend.key=element_blank(), legend.key.size=unit(2.5, "lines")) + 
      guides(colour=guide_legend(nrow=2, title=NULL))
    print(g3)
    dev.off()
    
    # Create Error rates plot: shows the values of four different error rates across the range of threshold probability values
    png(file=file.path(bccvl.env$outputdir, sprintf("%s-error-rates_%s.png", make.plot, species_algo_str)), width=480, height=480)
    g4 <- ggplot(errs[errs$measure %in% c("fpr", "fnr", "fdr", "fors"), ], 
                 aes(x=tpv, y=value, colour=measure, linetype=measure)) + 
      geom_line(size=1.2) + 
      ylim(0,1) +
      labs(title="Error rates plot\n", x="\nThreshold probability value", y="Error rate value\n") + 
      scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed"), labels=c("False Positive Rate  ", "False Negative Rate   ", "False Discovery Rate", "False Omission Rate")) +
      scale_colour_manual(values=c("#FAB334", "#D55E00", "#FAB334", "#D55E00"), labels=c("False Positive Rate  ", "False Negative Rate   ", "False Discovery Rate", "False Omission Rate")) + 
      theme(axis.text = element_text(family="Arial", size=rel(1.5)), axis.title = element_text(family="Arial", size=rel(1.5)), plot.title = element_text(family="Arial", size=rel(2)), legend.text = element_text(family="Arial", size=rel(1.5)), legend.position="top", legend.key=element_blank(), legend.key.size=unit(2.5, "lines")) +
      guides(colour=guide_legend(nrow=2, title=NULL), linetype=guide_legend(nrow=2, title=NULL))
    print(g4)
    dev.off()
    
    # Create ROC plot 
    png(file=file.path(bccvl.env$outputdir, sprintf("%s-ROC_%s.png", make.plot, species_algo_str)), width=480, height=480)
    xmax1 = min(round((max(temp$fpr) + 0.2)/0.1)*0.1, 1)
    xpos = max(xmax1/2, 0.1) 
    g5 <- ggplot(temp, aes(x=fpr, y=tpr)) + 
      geom_line(size=1.2) + 
      ylim(0,1) +
      xlim(0, xmax1) +
      geom_abline(intercept=0, slope=1, colour="grey") + 
      labs(x="\nFalse Positive Rate (1-Specificity)", y="True Positive Rate (Sensitivity)\n") +
      ggtitle(paste("ROC plot")) +
      annotate(geom = "text", x = xpos, y = 0.1, label = paste("AUC = ", auc.d$auc), size = 6) +
      theme(axis.text = element_text(family="Arial", size=rel(1.5)), axis.title = element_text(family="Arial", size=rel(1.5)), plot.title = element_text(family="Arial", size=rel(2)), legend.text = element_text(family="Arial", size=rel(1.5)))
    print(g5)
    dev.off()
    
    # Create evaluation stats table with values for optimum tpv value for L.diag (= maximize sum of TPR + TNR). 
    all.stats <- data.frame(list(tpv=list.tpv, L.diag=temp$L.diag, tpr=tpr, tnr=tnr, fpr=fpr, fnr=fnr, fdr=fdr, fors=fors, ppv=ppv, npv=npv, kappa=kappa, tss=tss, bs=bs, csi=csi, ets=ets, or=or, acc=acc, mcr=mcr))   
    all.stats <- round(all.stats, digits = 3)
    
    # Addition Number 2, 23 July 2018
    max.TPR.TNR <- all.stats[pick1min(all.stats$L.diag), ] # select row with all stats for maximum value of L.diag
    rownames(max.TPR.TNR) <- c("maximize sum TPR and TNR, minimize sum FPR and FNR")
    # End Addition Number 2, 23 July 2018

    # stats.table <- rbind(max.TPR.TNR, TPR.eq.TNR, max.Kappa) 
    stats.table <- max.TPR.TNR
    stats.table$L.diag <- NULL
    names(stats.table) <- c("Optimum threshold value:", "True Positive Rate (TPR)", "True Negative Rate (TNR)", "False Positive Rate (FPR)", "False Negative Rate (FNR)", 
                            "False Discovery Rate (FDR)", "False Omission Rate (FOR)", "Positive Predictive Value (PPV)", "Negative Predictive Value (NPV)", 
                            "Cohen's Kappa", "True Skill Statistic (TSS)", "Bias Score (BS)", "Critical Success Index (CSI)", "Equitable Threat Score (ETS)",
                            "Odds-Ratio (OR)","Accuracy", "Misclassification Rate")
    eval.stats <- t(stats.table) # transpose table

     # Create Loss function plot: shows the values of different loss functions across the range of threshold probability values
    png(file=file.path(bccvl.env$outputdir, sprintf("%s-loss-functions_%s.png", make.plot, species_algo_str)), width=480, height=480)
    g6 <- ggplot(errs[errs$measure %in% rev(c("L.diag", "L.pred", "L.all", "L.eq.diag")), ], 
                 aes(x=tpv, y=value, colour=measure)) + 
      geom_line(size=1.2) + 
      ylim(0,1) +
      labs(title="Loss function plot\n", x="\nThreshold probability value", y="Loss function value\n") +
      scale_colour_manual(values=c("#48D1CC", "#9F79EE", "#EE9572", "#FF3E96"), labels=c("Maximize TPR + TNR   ", "Maximize PPV + NPV   ", "Balance all errors", "TPR = TNR")) + 
      theme(axis.text = element_text(family="Arial", size=rel(1.5)), axis.title = element_text(family="Arial", size=rel(1.5)), plot.title = element_text(family="Arial", size=rel(2)), legend.text = element_text(family="Arial", size=rel(1.5)), legend.position="top", legend.key=element_blank(), legend.key.size=unit(2.5, "lines")) + 
      guides(colour=guide_legend(nrow=2, title = NULL))
    print(g6)
    dev.off()
    
    # Create Loss functions-intervals plot within 5% of the best value
    rangeperf$type.of.loss <- factor(rangeperf$type.of.loss, levels=(c("diag", "pred", "all", "eq.diag")))
    png(file=file.path(bccvl.env$outputdir, sprintf("%s-loss-intervals_%s.png", make.plot, species_algo_str)), width=480, height=480)
    g7 <- ggplot(rangeperf, aes(x=type.of.loss, y=best, ymin=lower, ymax=upper, colour=type.of.loss)) + 
      geom_pointrange(size=1.2) + 
      geom_line(size=1.2) +
      coord_flip() + 
      ylim(0,1) +
      scale_x_discrete(limits=c("diag","pred", "all", "eq.diag")) +
      scale_colour_manual(values=c("#48D1CC", "#9F79EE", "#EE9572", "#FF3E96"), labels=c("Maximize TPR + TNR   ", "Maximize PPV + NPV   ", "Balance all errors   ", "TPR = TNR")) + 
      labs(title="Range of threshold probability value \nwithin 5% of minimum per loss\n", x="Type of loss function\n", y="\nThreshold probability value") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(family="Arial", size=rel(1.5)), axis.title = element_text(family="Arial", size=rel(1.5)), plot.title = element_text(family="Arial", size=rel(2)), legend.text = element_text(family="Arial", size=rel(1.5)), legend.position="top", legend.key=element_blank(), legend.key.size=unit(2.5, "lines")) + 
      guides(colour=guide_legend(nrow=2, title = NULL))
    print(g7)
    dev.off()
  }
  
  return(list(performance=temp, stats=eval.stats, loss.summary=loss.table)) 
}

dev.save <- function(fileroot, ext=".pdf") {
  if (ext==".eps") {dev.copy2eps(file=paste(fileroot,ext,sep="."))} else {dev.copy2pdf(file=paste(fileroot,"pdf",sep="."))}
}

bccvl.createMarginalResponseCurves <- function(out.model, model.name, species_algo_str) {
  # Get the enviromental variables and values used to create the model
  if (model.name == "brt") {
    model.values = matrix(out.model$data$x, ncol=length(out.model$var.names))
    env.vars = out.model$var.names
  } else if (model.name %in% c("geoIDW", "voronoiHull")) {
    model.values = rbind(out.model@presence, out.model@absence)
    env.vars = colnames(model.values)
  } else {
    model.values = out.model@presence
    env.vars = colnames(model.values)
  }
  
  if (!(length(model.values)==0)) {
    # Create a matrix to hold average values for each environmental variable
    mean.values = matrix(data = NA, nrow = 100, ncol = length(env.vars))
    colnames(mean.values) = env.vars
    # For each variable, populate the column with the mean value
    for (i in 1:ncol(mean.values)) {
      mean.values[,i] = rep(mean(model.values[,i], na.rm=TRUE), 100)
    }

    # plot 18 response curves per page
    curvesPerPage = 6*3       # No of rows X No of columns
    for (i in 0:((ncol(mean.values)-1)/curvesPerPage)) {
        png(file=file.path(bccvl.env$outputdir, sprintf("response_curve_%s_p%d.png", species_algo_str, i)), width=700, height=900)
        par(mfrow = c(6,3)) # No of rows X No of columns

        # Allow each environmental variable to vary, keeping other variable values at average, and predict suitability
        rcurves = list()
        for (j in ((i*curvesPerPage + 1):min((i+1)*curvesPerPage, ncol(mean.values)))) {
            range.values = seq(min(model.values[,j], na.rm=TRUE), max(model.values[,j], na.rm=TRUE), length.out=100)
            temp.data = mean.values
            temp.data[,j] = range.values
            if (model.name == "brt") {
                colnames(temp.data) = env.vars
                new.predictions = predict(out.model, as.data.frame(temp.data), n.trees = out.model$gbm.call$best.trees, type="response")
            } else {
                new.predictions = predict(out.model, temp.data)
            }

            # create separate file for each response curve
            save.name = env.vars[j]
            plot(range.values, new.predictions, ylim=c(0,1), xlab="", ylab="", main=save.name, type="l")
            rug(model.values[,j])

            # save the response curve for later use
            df1 = data.frame(range.values, new.predictions)
            names(df1) <- c(save.name, "")
            rcurves[[save.name]] = df1
        }
        dev.off()

        # Save each response curve 
        for (k in 1:length(rcurves))
        {
          ename = env.vars[k + i*curvesPerPage]
          png(file=file.path(bccvl.env$outputdir, sprintf("%s_response_curve_%s.png", ename, species_algo_str)))
          plot(rcurves[[ename]], ylim=c(0,1), xlab="", ylab="", main=ename, type="l")
          rug(model.values[, k + i*curvesPerPage])
          dev.off()
        }
        rcurves = NULL
    }
  } else {
    write(paste(species_algo_str, ": Cannot create response curves from", model.name, "object", sep=" "), stdout())
  }
}

# function to calculate variable importance values for dismo models based on biomod2's correlation between predictions
# i.e., hold all but one predictor variable to its actual values, resample that one predictor and recalculate model predictions
bccvl.calculateVariableImpt <- function(out.model, model.name, num_samples, species_algo_str) {
  # EMG num_samples should be same as biomod.VarImport arg set in
  # 01.init.args.model.current.R
  
  # get the enviromental variables and values used to create the model
  # EMG this is duplicated from above, should be able to combine
  if (model.name == "brt") {
    model.values = matrix(out.model$data$x, ncol=length(out.model$var.names))
    env.vars = out.model$var.names
    colnames(model.values) = env.vars
  } else if (model.name %in% c("geoIDW", "voronoiHull")) {
    model.values = rbind(out.model@presence, out.model@absence)
    env.vars = colnames(model.values)
  } else {
    model.values = out.model@presence
    env.vars = colnames(model.values)
  }
  
  if (!(length(model.values)==0)) {
    # predict using actual values
    if (model.name == "brt") {
      actual.predictions = predict(out.model, as.data.frame(model.values), n.trees = out.model$gbm.call$best.trees, type="response")
    } else {
      actual.predictions = predict(out.model, model.values)
    }
    # create a table to hold the output
    varimpt.out = matrix(NA, nrow=length(env.vars), ncol=num_samples+2)
    dimnames(varimpt.out) = list(env.vars, c(paste("sample_", c(1:num_samples, "mean")), "percent"))
    # create a copy of the env data matrix
    sample.data = model.values
    # for each predictor variable
    for (p in 1:ncol(sample.data)) {
      # for each num_sample
      for (s in 1:num_samples) {
        # resample from that variables' values, keeping other variable values the same, and predict suitability
        sample.data[,p] = sample(x=sample.data[,p], replace=FALSE)
        # predict using sampled values
        if (model.name == "brt") {
          new.predictions = predict(out.model, as.data.frame(sample.data), n.trees = out.model$gbm.call$best.trees, type = "response")
        } else {
          new.predictions = predict(out.model, sample.data)
        }
        # calculate correlation between original predictions and new predictions
        varimpt.out[p,s] = 1-max(round(cor(x=actual.predictions, y=new.predictions, use="pairwise.complete.obs", method="pearson"), digits=3),0)
      }
    }
    # calculate mean variable importance, normalize to percentages, and write results
    varimpt.out[,num_samples+1] = round(rowMeans(varimpt.out, na.rm=TRUE), digits=3)
    varimpt.out[,num_samples+2] = round((varimpt.out[,num_samples+1]/sum(varimpt.out[,num_samples+1]))*100, digits=0)
    bccvl.write.csv(varimpt.out, name=sprintf("biomod2_like_VariableImportance_%s.csv", species_algo_str))
  } else {
    write(paste(species, ": Cannot calculate variable importance for ", model.name, "object", sep=" "), stdout())
  }
}

# function to calculate variable importance values for dismo models based on Maxent's decrease in AUC
# i.e., hold all but one predictor variable to its original values, resample that one predictor and recalculate model AUC
bccvl.calculatePermutationVarImpt <- function(out.model, model.eval,
                                              model.name, occur, bkgd, species_algo_str) {
  # get the enviromental variables and values used to create the model
  # EMG this is duplicated from above, should be able to combine or find an easier way to determine
  if (model.name == "brt") {
    model.values = matrix(out.model$data$x, ncol=length(out.model$var.names))
    env.vars = out.model$var.names
    colnames(model.values) = env.vars
  } else if (model.name %in% c("geoIDW", "voronoiHull")) {
    model.values = rbind(out.model@presence, out.model@absence)
    env.vars = colnames(model.values)
  } else {
    model.values = out.model@presence
    env.vars = colnames(model.values)
  }
  
  if (!(length(model.values)==0)) {
    # get the occurrence and background environmental data used to evaluate the model
    p.swd=occur
    a.swd=bkgd
    # get the AUC from the original model evaluation
    init.auc = round(model.eval@auc, digits=3)
    # create a table to hold the output
    permvarimpt.out = matrix(NA, nrow=length(env.vars), ncol=4)
    dimnames(permvarimpt.out) = list(env.vars, c("init.auc", "sample.auc", "change.auc", "percent"))
    permvarimpt.out[,"init.auc"] = rep(init.auc, length(env.vars))
    # create a copy of the occurrence and background environmental data
    sample.p = p.swd[,env.vars, drop=FALSE]
    sample.a = a.swd[,env.vars, drop=FALSE]
    # check for and remove any NA's present in the data
    no.na.sample.p = na.omit(sample.p);
    no.na.sample.a = na.omit(sample.a)
    if (nrow(no.na.sample.p) != nrow(sample.p)) {
      write(paste("bccvl.calculatePermutationVarImpt(): NA's were removed from presence data!"), stdout())
    }
    if (nrow(no.na.sample.a) != nrow(sample.a)) {
      write(paste("bccvl.calculatePermutationVarImpt(): NA's were removed from absence data!"), stdout())
    }
    # for each predictor variable
    for (v in 1:length(env.vars)) {
      # resample from that variables' values, keeping other variable values the same
      no.na.sample.p[,v] = sample(x=no.na.sample.p[,v], replace=FALSE)
      no.na.sample.a[,v] = sample(x=no.na.sample.a[,v], replace=FALSE)
      # re-evaluate model with sampled env values
      if (model.name == "brt") {
        sample.eval = dismo::evaluate(p=no.na.sample.p, a=no.na.sample.a, model=out.model, n.trees=out.model$gbm.call$best.trees, type="response")
      } else {
        sample.eval = dismo::evaluate(p=no.na.sample.p, a=no.na.sample.a, model=out.model)
      }
      # get the new auc
      permvarimpt.out[v,"sample.auc"] = round(sample.eval@auc, digits=3)
    }
    # calculate the difference in auc, normalize to percentages, and write results
    permvarimpt.out[,"change.auc"] = permvarimpt.out[,"init.auc"] - permvarimpt.out[,"sample.auc"]
    for (r in 1:nrow(permvarimpt.out)) {
      if (permvarimpt.out[r,"change.auc"] < 0) {  # EMG what if AUC increases?
        permvarimpt.out[r,"change.auc"] = 0
      }
    }
    permvarimpt.out[,"percent"] = round((permvarimpt.out[,"change.auc"]/sum(permvarimpt.out[,"change.auc"]))*100, digits=0)
    bccvl.write.csv(permvarimpt.out, name=sprintf("maxent_like_VariableImportance_%s.csv", species_algo_str))
  } else {
    write(paste(species, ": Cannot calculate maxent-like variable importance for ", model.name, "object", sep=" "), stdout())
  }
}

#########################################################################
### 4: Run the evaluation and save outputs
#########################################################################

### Evaluation of 'dismo' models and save outputs

bccvl.saveDISMOModelEvaluation <- function(model.name, model.obj, occur, bkgd, species.name) {
  species_algo_str = paste(species.name, model.name, sep="_")

  # evaluate model using dismo's evaluate
  if (model.name == "brt") {
    model.eval = dismo::evaluate(p=occur, a=bkgd, model=model.obj, n.trees=model.obj$gbm.call$best.trees, type="response")
  } else {
    model.eval = dismo::evaluate(p=occur, a=bkgd, model=model.obj)
  }
  # need predictions and observed values to create confusion matrices for accuracy statistics
  model.fit = c(model.eval@presence, model.eval@absence)
  model.obs = c(rep(1, length(model.eval@presence)), rep(0, length(model.eval@absence)))
  
  # Call the evaluation script
  res = performance.2D(model.obs, model.fit, species_algo_str, make.plot="dismo", kill.plot=F)
  bccvl.saveModelEvaluation(res$performance, res$stats, res$loss.summary, species_algo_str)
  
  # Create response curves
  bccvl.createMarginalResponseCurves(model.obj, model.name, species_algo_str)
  
  # Calculate variable importance (like biomod2, using correlations between predictions)
  bccvl.calculateVariableImpt(model.obj, model.name, 3, species_algo_str)
  
  # Calculate variable importance (like maxent, using decrease in AUC)
  bccvl.calculatePermutationVarImpt(model.obj, model.eval, model.name, occur, bkgd, species_algo_str)
  
  # Create HTML file with accuracy measures
  # bccvl.generateHTML()
}

### Evaluation of 'biomod2' models and save outputs

bccvl.saveBIOMODModelEvaluation <- function(loaded.names, biomod.model, species_algo_str) {

  evaluation = get_evaluations(biomod.model)

  # Get the model predictions and observed values. Predictions is a 4-dimensional array (Predictions, Algorithm, Model run, PseudoAbsence Run)
  predictions = get_predictions(biomod.model)
  total_models = length(dimnames(predictions)[[3]])
  
  obs = get_formal_data(biomod.model, "resp.var")
  # In case of pseudo-absences we might have NA values in obs so replace them with 0
  obs = replace(obs, is.na(obs), 0)
  
  for ( i in 1:total_models )
  {
    model_name = dimnames(predictions)[[3]][i]  # will be FULL or RUN1 for eg
    model_predictions = predictions[,,i,]
    
    if (sum(is.na(model_predictions)) == length(model_predictions)) 
    {
      # Warn that model n is being ignored. It most probably failed to build.
      warning(sprintf("Warning: Model %i failed to generate. Not generating stats", i), immediate.=T)
      next
    }

    res = performance.2D(obs, model_predictions / 1000, species_algo_str, make.plot=model_name, kill.plot=F)
    bccvl.saveModelEvaluation(res$performance, res$stats, res$loss.summary, species_algo_str)
    
    # get and save the variable importance estimates
    variableImpt = get_variables_importance(biomod.model)
    if (!is.na(variableImpt)) {
      #EMG Note this will throw a warning message if variables (array) are returned
      bccvl.write.csv(variableImpt, 
                      name=paste("variableImportance", model_name, species_algo_str, "csv", sep="."))
    } else {
      message("VarImport argument not specified during model creation!")
      #EMG must create the model with the arg "VarImport" != 0
    }
  }
  
  # save response curves (Elith et al 2005)
  for(name in loaded.names)
  {
    env_data = get_formal_data(biomod.model,"expl.var")
    png(file=file.path(bccvl.env$outputdir, sprintf("mean_response_curves_%s.png", name)))
    test <- response.plot2(models = name,
                           Data = env_data,
                           show.variables = get_formal_data(biomod.model,"expl.var.names"),
                           fixed.var.metric = "mean")
    dev.off()

    # save individual response curves
    for (envname in names(test))
    {
      png(file=file.path(bccvl.env$outputdir, sprintf("%s_mean_response_curves_%s.png", envname, name)))
      plot(test[[envname]], type='l', ylim=c(0, 1.0), main=envname, xlab="", ylab="")
      rug(env_data[[envname]])
      dev.off()
    }
  }
}

#
# ------------------------------------------
#

bccvl.savePdf <- function(..., filename, aspdf, outputdir=bccvl.env$outputdir)
{
  library("gridExtra")
  if (aspdf) 
  { 
    png(file=file.path(outputdir, paste(filename, 'png', sep=".")))
    grid.arrange(...)
    dev.off()  
  }
  else {     
    grid.arrange(...)
  }
}


# A special R function for variable importance plots based on biomod2 fitted model outcomes

# TODO: data1 passed in here may encode pseudo absence as NA ... we should convert that to 0
# TODO: could use checks for covars != NA in data1
bccvl.VIPplot <- function(fittedmodel=NULL, 
                         method=c("glm",,"cta","gam","ann", "rf", "gbm", "mars", "maxent"),  
                         cor.method=c("pearson","spearman"),
                         pdf=TRUE, biom_vi=FALSE,output.table=FALSE, data1, this.dir, filename)
{
  library("ggplot2")
  library("reshape2")
  library("mgcv")
  library("rpart")
  library("caret")
  library("ggdendro")

  # README notes:
  # (1) fittedmodel: the fitted model object obtained running the biomod2 function'BIOMOD_Modeling'.
  #     method is one of "glm","rpart","gam","ann", "rf", "gbm", "mars", "maxent.p"
  #     cor.method is one of "pearson","spearman"
  # (2) pdf=FALSE: the default setting - no pdf file generated;  pdf=TRUE: the VIP generated is
  #     saved as a pdf file in the working directory.
  # (3) biom_vi=FALSE: a function/algorithm other than the biomod2 inbuilt function 'variables_importance' 
  #     will be applied for evaluating/ranking the variable importance.
  # (4) output.table=FALSE: a .csv file which contains a table to display the glm model parameter
  #     estimates and the 95% confidence bounds for both raw data and scaled data will be generated
  #     if output.table=TRUE.
  # (5) cor.method=c("pearson","spearman"): the default "pearson" method measures only the linear 
  #     association among variables; the "spearman" method is a rank-based algorithm which is more
  #     robust measure for association among variables, e.g., the non-linear association will also be detected.
  # (6) Note that 'data1' is a dataframe with the response variable in the first column following
  #     by the predictor variables / enviornmental variables in the other columns;
  #     'data1' is needed for generating the VIP by the biomod2 inbuilt function
  #     'variables_importance(fittedmodel, data1)' and for calculate the AIC scores.
  #     Warning: Except for the glm algorithm, the 'data1' should include all predictor variables to make
  #     the variable importance ranking outcomes meaningful.
  # (7) this.dir specifies the route to access the biomod2 model (but not including the model name)
  # (8) filename to be saved without the file extension.
  # 

  data1$y.data[is.na(data1$y.data)] <- 0
  
 # extract the root of filenames used by biomod to save model results
 filenames <- dir(this.dir)
 #loc <- regexpr("_RUN[[:digit:]]", filenames[1])
 #fileroot <- substr(filenames[1], 1, loc-1)

# select the full model generated
 filekeep <-  paste(this.dir, "/", filenames[1], sep="")
 
 if (!is.na(match("glm",method))) 
 {
   working <- load(filekeep)
   fittedmodel <- get_formal_model(eval(parse(text=working)))

   fitteddata = fittedmodel$data  # these are the data used by biomod2 for model fitting
   nd = dim(fitteddata)[2]
   sub.data = fitteddata[,2:nd, drop=FALSE]
   # Can only scale numeric data
   cat.data = Filter(is.factor, sub.data)
   sub.data = Filter(is.numeric, sub.data)
   RespV = fitteddata[,1, drop=FALSE]
   rescaled.data <- scale(sub.data)

   # attributes(rescaled.data)
   all.data <- cbind(RespV, as.data.frame(rescaled.data), cat.data)

   # head(all.data); dim(all.data)
   rescaled.glm <- glm(fittedmodel$formula, data=all.data, family=binomial)

   scaled = as.data.frame(cbind(coef(rescaled.glm), confint(rescaled.glm)))
   scaled = scaled[-1,]
   raw = as.data.frame(cbind(coef(fittedmodel), confint(fittedmodel)))
   raw = raw[-1,]

   #  variable importance plot in terms of relative effect size
   #  the relative effect size is defined as the coefficient estimated based on scaled data
   nx = length(coef(fittedmodel)[-1])
   df1 = as.data.frame(cbind(1:nx,round(scaled,3)))

   names(df1) = c("xvar", "meanest", "lower", "upper")
   df1$xvar = factor(df1$xvar, labels = rownames(df1))

   p1 <- ggplot(df1, aes(x=xvar, y=meanest)) + geom_hline(yintercept = 0) + labs(x=" ")

   ps = p1 + geom_errorbar(aes(ymin=lower,ymax=upper),lwd=0.8,width=0.25) + 
      labs(y="relative effect size") +  labs(title="         scaled data") + coord_flip()

   df2 = as.data.frame(cbind(1:nx,round(raw,3)))
   names(df2) = c("xvar", "meanest", "lower", "upper")
   df2$xvar = factor(df2$xvar, labels = rownames(df2))

   if (output.table) 
   {
      df1t = df1; df2t = df2
      names(df1t) = c("x.var", "coeff.est.scaled", "lower", "upper")
      names(df2t) = c("x.var", "coeff.est.raw", "lower", "upper")
      dfout = cbind(df2t,df1t)
      write.csv(dfout,file=paste(filekeep,"paraest_out.csv",sep="_"),row.names=FALSE)
   }

   #  the heatmap in terms of correlation among numerical predictor variables
   rescdata = Filter(is.numeric, rescaled.glm$data[,-1, drop=FALSE])

   if("spearman" %in% cor.method) {
       xx = cor(rescdata, method="spearman")
   } else if("pearson" %in% cor.method)  {
       xx = cor(rescdata)
   }

   lower_tri <- xx
   lower_tri[upper.tri(lower_tri)] <- NA
  
   xx.ml <- melt(lower_tri,na.rm=TRUE)  #the argument 'na.rm=TRUE' seems not working)

   corx = xx.ml[,3]
   rm = which(is.na(corx)==TRUE)
   xx.ml = xx.ml[-rm,]

   pheat <- ggplot(xx.ml, aes(X1, X2)) + geom_tile(aes(fill = value), colour="black") + 
      scale_fill_gradient2(low = "green4", high = "violetred", mid="white", 
      midpoint=0, limit=c(-1,1)) + labs(y=" ") + theme_minimal() +
      scale_x_discrete(limits=rownames(xx)) + scale_y_discrete(limits=colnames(xx)) + coord_fixed() +
      theme(axis.title.x=element_blank(),legend.position = "bottom", axis.text.x=element_text(angle=-90)) +
      guides(fill=guide_legend(title="correlation"))

   # Save as variable correlation plot.
   filename1 = sub("vip_plot", "variable_correlations", filename)
   bccvl.savePdf(pheat, ncol=1, nrow=1, filename=filename1, aspdf=pdf)

   # variable importance plot in terms of AIC scores which represent the information loss,
   # e.g., the AIC score of a predictor variable representing the information loss 
   # if this variable is not included in the selected model.

   nd = dim(data1)[2]

   RespV1 = data1[,1]; subdata1 = data1[,2:nd, drop=FALSE]
   glm.all = glm(formula = RespV1 ~ ., family = binomial, data = subdata1)

   Xaic = NULL
   for (i in 1:(nd-1))
   {
      subdf = subdata1[,-i, drop=FALSE]
      glm.one = glm(formula = RespV1 ~ . , family = binomial, data = subdf)
      Xaic = c(Xaic,AIC(glm.one)) 
   }

   relaAIC = round(Xaic - AIC(glm.all),2)  
   nx = length(relaAIC)
   dfa = as.data.frame(cbind(1:nx,relaAIC))
   dfa$V1 = factor(dfa$V1, labels = rev(names(subdata1)))
   pa <- ggplot(dfa, aes(x=V1, y=rev(relaAIC))) + labs(x="predictor variables") + 
   labs(y="AIC score for information loss") + labs(title="AIC approach")
  
   ppa = pa + geom_col(alpha=0.6,col="blue") + coord_flip()

   # the variable importance plot using the inbuilt biomod2 function 'variables_importance'
   vi_biomod = variables_importance(fittedmodel,data=subdata1)$mat
   nx = length(vi_biomod)
   dfvi = as.data.frame(cbind(1:nx,vi_biomod[,1]))
   dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(vi_biomod)))
   pv <- ggplot(dfvi, aes(x=V1, y=rev(vi_biomod[,1]))) + labs(x="predictor variables") + 
      labs(y="variable importance score") + labs(title="biomod2 function 'variables_importance'")
 
   ppv = pv + geom_col(alpha=0.6,col="green4") + coord_flip()

   # Save as variable relative contribution plot.
   filename1 = sub("vip_plot", "variable_relative_contribution", filename)
   if (biom_vi) {
     bccvl.savePdf(ps, ppv, ncol=2, nrow=1, filename=filename1, aspdf=pdf)
   }
   else {
     bccvl.savePdf(ps, ppa, ncol=2, nrow=1, filename=filename1, aspdf=pdf)
   }
 }

 if (!is.na(match("cta",method))) 
 {
   working <- load(filekeep)
   fittedmodel <- get_formal_model(eval(parse(text=working)))

# variable importance is part of the model fitting outcomes with 'rpart' algorithm and
# this information can be used for generating the variable importance plot

   varimp0 = fittedmodel$variable.importance
   nx = length(varimp0)
   df0 = as.data.frame(cbind(1:nx,varimp0))
   df0$V1 = factor(df0$V1, labels = rev(names(varimp0)))

   p <- ggplot(df0, aes(x=V1, y=rev(varimp0))) + labs(x=" ") + 
      labs(y="variable importance score") + labs(title="part of the 'rpart' model output")

   pp0 = p + geom_col(alpha=0.6,col="blue") + coord_flip()

   ddata <- dendro_data(fittedmodel)
   ppt = ggplot() + 
      geom_segment(data = ddata$segments, aes(x = x, y = y, xend = xend, yend = yend)) + 
      geom_text(data = ddata$labels, aes(x = x, y = y, label = label), size = 3, vjust = 0) +
      geom_text(data = ddata$leaf_labels, aes(x = x, y = y, label = label), size = 3, vjust = 1) +
      theme_dendro()

   # variable importance plot using the inbuilt biomod2 function 'variables_importance'

   nd = dim(data1)[2]
   subdata1 = data1[,2:nd, drop=FALSE]
   vi_biomod = variables_importance(fittedmodel,data=subdata1)$mat
   nx = length(vi_biomod)
   dfvi = as.data.frame(cbind(1:nx,vi_biomod[,1]))
   dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(vi_biomod)))
   pv <- ggplot(dfvi, aes(x=V1, y=rev(vi_biomod[,1]))) + labs(x="predictor variables") + 
      labs(y="variable importance score") + labs(title="biomod2 function 'variables_importance'")
   ppv = pv + geom_col(alpha=0.6,col="green4") + coord_flip()

   if (biom_vi) 
   {
     bccvl.savePdf(ppt, ppv, ncol=2, nrow=1, filename=filename, aspdf=pdf)
   }
   else
   {
     bccvl.savePdf(ppt, pp0, ncol=2, nrow=1, filename=filename, aspdf=pdf)
   }
 }


 if (!is.na(match("gam",method))) 
 {
   working <- load(filekeep)
   fittedmodel <- get_formal_model(eval(parse(text=working)))

   if (biom_vi) 
   {
     # variable importance plot using the inbuilt biomod2 function 'variables_importance'
     nd = dim(data1)[2]  
     RespV1 = data1[,1]; subdata1 = data1[,2:nd, drop=FALSE]
   
     vi_biomod = variables_importance(fittedmodel,data=subdata1)$mat
     nx = length(vi_biomod)
     dfvi = as.data.frame(cbind(1:nx,vi_biomod[,1]))
     dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(vi_biomod)))
     pv <- ggplot(dfvi, aes(x=V1, y=rev(vi_biomod[,1]))) + labs(x="predictor variables") + 
        labs(y="variable importance score") + labs(title="biomod2 function 'variables_importance'")
     ppv = pv + geom_col(alpha=0.6,col="green4") + coord_flip()
    
     bccvl.savePdf(ppV, ncol=1, nrow=1, filename=filename, aspdf=pdf)
   }  # end of 'if(biom_vi=TRUE)' 
   else 
   {
     # variable importance plot following the AIC approach
     nd = dim(data1)[2]  
     RespV1 = data1[,1]
     subdata1 = data1[,2:nd, drop=FALSE]

     # gam function cannot take categorical data, so exclude categorical data.
     subdata1 = Filter(is.numeric, subdata1)

     xname = names(subdata1)
     sname = paste("s(", xname, ")",sep="")

     gamformu.all <- as.formula(paste("RespV1 ~ 1 +", paste(sname, collapse= "+")))
     gam.all = gam(formula = gamformu.all, family = binomial, data = subdata1)

     Xaic = NULL
     nd = dim(subdata1)[2]
     for (i in 1:nd)
     {
        subdf = subdata1[, -i, drop=FALSE]
        xname1 = names(subdf)
        sname1 = paste("s(", xname1, ")",sep="")
        gamformu1 <- as.formula(paste("RespV1 ~ 1 +", paste(sname1, collapse= "+")))
        gam.one = gam(formula = gamformu1, family = binomial, data = subdf)
        Xaic = c(Xaic,AIC(gam.one)) 
     }

     relaAIC = round(Xaic - AIC(gam.all),2)
     nx = length(relaAIC)
     dfa = as.data.frame(cbind(1:nx,relaAIC))
     dfa$V1 = factor(dfa$V1, labels = rev(names(subdata1)))
     pa <- ggplot(dfa, aes(x=V1, y=rev(relaAIC))) + labs(x="predictor variables") + 
        labs(y="AIC score for information loss") + labs(title="AIC approach")
    
     ppa = pa + geom_col(alpha=0.6,col="blue") + coord_flip()

     bccvl.savePdf(ppa, ncol=1, nrow=1, filename=filename, aspdf=pdf)

   }
 }

 if (!is.na(match("ann",method))) 
 {
   working <- load(filekeep)
   fittedmodel <- get_formal_model(eval(parse(text=working)))

   # variable importance plot using the inbuilt biomod2 function 'variables_importance'
   nd = dim(data1)[2]  
   RespV1 = data1[,1]; subdata1 = data1[,2:nd]
 
   vi_biomod = variables_importance(fittedmodel,data=subdata1)$mat
   nx = length(vi_biomod)
   dfvi = as.data.frame(cbind(1:nx,vi_biomod[,1]))
   dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(vi_biomod)))
   pv <- ggplot(dfvi, aes(x=V1, y=rev(vi_biomod[,1]))) + labs(x="predictor variables") + 
      labs(y="variable importance score") + labs(title="biomod2 function 'variables_importance'")
   ppv = pv + geom_col(alpha=0.6,col="green4") + coord_flip()
  
   bccvl.savePdf(ppv, ncol=1, nrow=1, filename=filename, aspdf=pdf)
 }

 if (!is.na(match("mars",method))) 
 {
   # Note that the model class generated from MARS algorithm is not supported by 
   #  the inbuilt biomod2 function 'variables_importance', neither the AIC approach is applicable.
   # However, the function 'varImp' in package 'caret' accept the MARS model object for estimating
   #  the variable importance. GCV = generalized cross validation
   working <- load(filekeep)
   fittedmodel <- get_formal_model(eval(parse(text=working)))

   # variable importance plot using the inbuilt function 'varImp' from package 'caret'
   nd = dim(data1)[2]  
   RespV1 = data1[,1]; subdata1 = data1[,2:nd]
 
   var_imp = varImp(fittedmodel)
   nx = length(var_imp[,1])
   dfvi = as.data.frame(cbind(1:nx,var_imp[,1]))
   dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(var_imp)))
   pv <- ggplot(dfvi, aes(x=V1, y=rev(var_imp[,1]))) + labs(x="predictor variables") + 
      labs(y="relative reduction in GCV") + labs(title="function 'varImp' in package 'caret'")
   ppv = pv + geom_col(alpha=0.6,col="red") + coord_flip()
  
   bccvl.savePdf(ppv, ncol=1, nrow=1, filename=filename, aspdf=pdf)
 }
   
  
 if (!is.na(match("gbm",method))) 
 {
   working <- load(filekeep)
   fittedmodel <- get_formal_model(eval(parse(text=working)))
 
   # variable importance plot using the inbuilt biomod2 function 'variables_importance'
   nd = dim(data1)[2]  
   RespV1 = data1[,1]; subdata1 = data1[,2:nd]
 
   vi_biomod = variables_importance(fittedmodel,data=subdata1)$mat
   nx = length(vi_biomod)
   dfvi = as.data.frame(cbind(1:nx,vi_biomod[,1]))
   dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(vi_biomod)))
   pv <- ggplot(dfvi, aes(x=V1, y=rev(vi_biomod[,1]))) + labs(x="predictor variables") + 
      labs(y="variable importance score") + labs(title="biomod2 function 'variables_importance'")
   ppv = pv + geom_col(alpha=0.6,col="green4") + coord_flip()
  
   bccvl.savePdf(ppv, ncol=1, nrow=1, filename=filename, aspdf=pdf)
 }


 if (!is.na(match("rf",method))) 
 {
   working <- load(filekeep)
   fittedmodel <- get_formal_model(eval(parse(text=working)))
 
   # Random forests (rf) provide an improvement over bagged trees by way of a small tweak
   #  that decorrelates the trees.
   # Note that the variable importance plot using the inbuilt biomod2 function 'variables_importance'
   #  does not seem working with Random forests algorithm.  On the other hand, however, the fitted
   #  rf model object contains the variable importance information which is measured by the mean derease
   #  in Gini index (expressed relative to the maximum).  
   # While RSS is used for measuring the regression tree model performance, the Gini index is used for 
   #  measuring the classification tree model performance and Gini index is a measure of total variance
   #  across the K classes.

   nd = dim(data1)[2]  
   RespV1 = data1[,1]; subdata1 = data1[,2:nd]
 
   out.rf = fittedmodel$importance
   rfImp = out.rf[,1]
   nx = length(rfImp)

   dfrf = as.data.frame(cbind(1:nx,rfImp))
   dfrf$V1 = factor(dfrf$V1, labels = rev(names(rfImp)))
   prf <- ggplot(dfrf, aes(x=V1, y=rev(rfImp))) + labs(x="predictor variables") + 
       labs(y="mean decrease in Gini index") + labs(title="part of rf model fitting outputs")
   pprf = prf + geom_col(alpha=0.6,col="blue") + coord_flip()
  
   bccvl.savePdf(pprf, ncol=1, nrow=1, filename=filename, aspdf=pdf)                      
 }  
 

 if (!is.na(match("maxent",method))) 
 {
   # variable importance is part of the model fitting outcomes with 'MAXENT.Phillips' algorithm
   if (regexpr("_outputs", filekeep) < 0)
   {
     working <- paste(filekeep,"_outputs/maxentResults.csv",sep="")
   }
   else
   {
     working <- paste(filekeep,"/maxentResults.csv",sep="")
   }

   df.P = read.csv(working)

   the.data <- data1[,-1]
   the.data <- Filter(is.numeric, the.data)
   nx = dim(the.data)[2]    #decide the number of the predictor variables
   yp = as.numeric(df.P[,8:(7+nx)])

   dfp = as.data.frame(cbind(1:nx,yp))
   dfp$V1 = factor(dfp$V1, labels = rev(names(the.data)))
   p <- ggplot(dfp, aes(x=V1, y=rev(yp))) + labs(x="predictor variables") + 
      labs(y="variable relative contribution (%)") + labs(title="maxent algorithm")
   pp = p + geom_col(alpha=0.6,col="red") + coord_flip()

   # Save as variable relative contribution plot.
   filename1 = sub("vip_plot", "variable_relative_contribution", filename)
   bccvl.savePdf(pp, ncol=1, nrow=1, filename=filename1, aspdf=pdf)

   #  the heatmap in terms of correlation among predictor variables
   if(cor.method=="pearson")  xx = cor(the.data)
   if(cor.method=="spearman") xx = cor(the.data,method="spearman")
 
    get_lower_tri<-function(cormat)
    {
      cormat[upper.tri(cormat)] <- NA
      return(cormat) 
    }

   lower_tri = get_lower_tri(xx)
  
   xx.ml <- melt(lower_tri,na.rm=TRUE)  #the argument 'na.rm=TRUE' seems not working)

   corx = xx.ml[,3]
   rm = which(is.na(corx)==TRUE)
   xx.ml = xx.ml[-rm,]

   pheat <- ggplot(xx.ml, aes(X1, X2)) + geom_tile(aes(fill = value), colour="black") + 
      scale_fill_gradient2(low = "green4", high = "violetred", mid="white", midpoint=0, limit=c(-1,1)) +
      scale_x_discrete(limits=rownames(xx)) + scale_y_discrete(limits=colnames(xx)) + coord_fixed() +
      labs(y=" ") + theme_minimal() +
      theme(axis.title.x=element_blank(),legend.position = "bottom", axis.text.x=element_text(angle=-90)) +
      guides(fill=guide_legend(title="correlation"))
   # Save as variable correlation plot.
   filename1 = sub("vip_plot", "variable_correlations", filename)
   bccvl.savePdf(pheat, ncol=1, nrow=1, filename=filename1, aspdf=pdf)
 }
}


########################################################
# Patch Biomod2 plot projection function
#  taken from: biomod2-3.1-64 : https://github.com/cran/biomod2/blob/3.1-64/R/BiomodClass.R#L1586
########################################################
require('rasterVis')
setMethod('plot', signature(x='BIOMOD.projection.out', y="missing"),
          function(x,col=NULL, str.grep=NULL, on_0_1000=TRUE){
            models_selected <- x@models.projected 
            if(length(str.grep)){
              models_selected <- grep(paste(str.grep,collapse="|"), models_selected,value=T)
            } 
            
            if(!length(models_selected)) stop("invalid str.grep arg")

            if(class(x@proj) == "BIOMOD.stored.raster.stack"){
              require(rasterVis)

              my.scale = ifelse(on_0_1000, 1, 1000)
              ## define the breaks of the color key
              my.at <- seq(0,1000,by=100) / my.scale
              ## the labels will be placed vertically centered
              my.labs.at <- seq(0,1000,by=250) / my.scale
              ## define the labels
              my.lab <- seq(0,1000,by=250) / my.scale
              ## define colors
              #               my.col <- colorRampPalette(c("red4","orange4","yellow4","green4"))(100)
              my.col <- colorRampPalette(c("grey90","yellow4","green4"))(100)

              ## try to use levelplot function
              try_plot <- try(
                levelplot(get_predictions(x, full.name=models_selected),
                          at=my.at, margin=T, col.regions=my.col,
                          main=paste(x@sp.name,x@proj.names,"projections"),
                          colorkey=list(labels=list(
                            labels=my.lab,
                            at=my.labs.at)))
              ) 
              if(! inherits(try_plot,"try-error")){ ## produce plot
                print(try_plot)
              } else{## try classical plot
                cat("\nrasterVis' levelplot() function failed. Try to call standard raster plotting function.",
                    "It can lead to unooptimal representations.",
                    "You should try to do it by yourself extracting predicions (see : get_predictions() function)", fill=options()$width)
                try_plot <- try(
                  plot(get_predictions(x, full.name=models_selected))
                )
              } 

              if(inherits(try_plot,"try-error")){ # try classical plot
                cat("\n Plotting function failed.. You should try to do it by yourself!")
              }
              
            } else if(class(x@proj) == "BIOMOD.stored.array"){
              if(ncol(x@xy.coord) != 2){
                cat("\n ! Impossible to plot projections because xy coordinates are not available !")
              } else {
                multiple.plot(Data = get_predictions(x, full.name=models_selected, as.data.frame=T), coor = x@xy.coord)
              }
            } else {cat("\n !  Biomod Projection plotting issue !", fill=.Options$width)}
          })

####
##
##  INPUT:
##
##  occur.data ... filename for occurence data
##  absen.data  ... filename for absence data
##  enviro.data.current ... list of filenames for climate data
##  enviro.data.type    ... continuous
##  opt.tails ... predict parameter
##
##  outputdir ... root folder for output data

#define the working directory
#scriptdir = normalizePath(bccvl.env$scriptdir)
#inputdir =  normalizePath(bccvl.env$inputdir)
#outputdir =  normalizePath(bccvl.env$outputdir)


# extract params
# define the lon/lat of the observation records -- 2 column matrix of longitude and latitude
occur.data = bccvl.params$species_occurrence_dataset$filename
occur.species = bccvl.params$species_occurrence_dataset$species
month.filter = bccvl.params$species_filter
#define the the lon/lat of the background / psuedo absence points to use -- 2 column matrix of longitude and latitude
absen.data = bccvl.params$species_absence_dataset$filename
#define the current enviro data to use
# rename file if filename contains 'asc' due to bug in Biomod2 function.
# TODO: Don't need to rename when Biomod2 bug is fixed.
enviro.data.current = lapply(bccvl.params$environmental_datasets,
                             function(x) {
                                fname = x$filename
                                if (file_ext(fname) != 'asc') {
                                   fname = file.path(dirname(fname), sub('asc', 'asd', basename(fname)))
                                   file.rename(x$filename, fname)
                                }
                                return(fname)
                             }
)
#type in terms of continuous or categorical
enviro.data.type = lapply(bccvl.params$environmental_datasets, function(x) x$type)
#layer names for the current environmental layers used
enviro.data.layer = lapply(bccvl.params$environmental_datasets, function(x) x$layer)
#geographic constraints.
enviro.data.constraints = readLines(bccvl.params$modelling_region$filename)
#Indicate to generate and apply convex-hull polygon of occurrence dataset to constraint
enviro.data.generateCHall = ifelse(is.null(bccvl.params$generate_convexhull), FALSE, as.logical(bccvl.params$generate_convexhull))
#Indicate whether to generate unconstraint map or not. True by default
enviro.data.genUnconstraintMap = ifelse(is.null(bccvl.params$unconstraint_map), TRUE, as.logical(bccvl.params$unconstraint_map))
# resampling (up / down scaling) if scale_down is TRUE, return 'lowest'
enviro.data.resampling = ifelse(is.null(bccvl.params$scale_down) ||
                                as.logical(bccvl.params$scale_down),
                                'highest', 'lowest')

############### BIOMOD2 Models ###############
#
# general parameters to perform any biomod modelling
#
biomod.NbRunEval = bccvl.params$nb_run_eval  # default 10; n-fold cross-validation; ignored if DataSplitTable is filled
biomod.DataSplit = bccvl.params$data_split # default 100; % for calibrating/training, remainder for testing; ignored if DataSplitTable is filled
biomod.Yweights = NULL #response points weights
biomod.Prevalence = bccvl.params$prevalence #either NULL (default) or a 0-1 numeric used to build "weighted response weights"
biomod.VarImport = bccvl.params$var_import # default 0; number of resampling of each explanatory variable to measure the relative importance of each variable for each selected model
#EMG this parameter needs to be specified in order to get VariableImportance metrics during model evaluation
biomod.models.eval.meth = c("KAPPA", "TSS", "ROC" ,"FAR", "SR", "ACCURACY", "BIAS", "POD", "CSI", "ETS") #vector of evaluation metrics
biomod.rescal.all.models = bccvl.params$rescale_all_models #if true, all model prediction will be scaled with a binomial GLM
biomod.do.full.models = bccvl.params$do_full_models #if true, models calibrated and evaluated with the whole dataset are done; ignored if DataSplitTable is filled
biomod.modeling.id = bccvl.params$modeling_id  #character, the ID (=name) of modeling procedure. A random number by default
# biomod.DataSplitTable = NULL #a matrix, data.frame or a 3D array filled with TRUE/FALSE to specify which part of data must be used for models calibration (TRUE) and for models validation (FALSE). Each column correspund to a "RUN". If filled, args NbRunEval, DataSplit and do.full.models will be ignored
# EMG Need to test whether a NULL values counts as an argument
biomod.species.name = occur.species # used for various path and file name generation
projection.name = "current"  #basename(enviro.data.current)
species_algo_str = ifelse(is.null(bccvl.params$subset),
                          sprintf("%s_maxent", occur.species),
                          sprintf("%s_maxent_%s", occur.species, bccvl.params$subset))


# model-specific arguments to create a biomod model
model.options.maxent <- list(
  path_to_maxent.jar = Sys.getenv("MAXENT"), #The path to maxent.jar file (the working directory by default)
  #memory_allocated = bccvl.params$memory_allocated #The amount of memory (in Mo) reserved for java to run MAXENT. should be 64, 128, 256, 512, 1024, 2048... or NULL if you want to use default java memory limitation parameter.
  maximumiterations = bccvl.params$maximumiterations,
  visible = FALSE,
  linear = bccvl.params$linear,
  quadratic = bccvl.params$quadratic,
  product = bccvl.params$product,
  threshold = bccvl.params$threshold,
  hinge = bccvl.params$hinge,
  lq2lqptthreshold = bccvl.params$lq2lqptthreshold,
  l2lqthreshold = bccvl.params$l2lqthreshold,
  hingethreshold = bccvl.params$hingethreshold,
  beta_threshold = bccvl.params$beta_threshold,
  beta_categorical = bccvl.params$beta_categorical,
  beta_lqp = bccvl.params$beta_lqp,
  beta_hinge = bccvl.params$beta_hinge,
  betamultiplier = bccvl.params$betamultiplier,
  defaultprevalence = bccvl.params$defaultprevalence
)

############### BIOMOD2 Models ###############
#
# general parameters to project any biomod modelling
#
#modeling.output #"BIOMOD.models.out" object produced by a BIOMOD_Modeling run
#new.env #a set of explanatory variables onto which models will be projected; must match variable names used to build the models
#proj.name #a character defining the projection name (a new folder will be created with this name)
# pseudo absences
biomod.PA.nb.rep = 0
biomod.PA.nb.absences = 0

biomod.xy.new.env = NULL #optional coordinates of new.env data. Ignored if new.env is a rasterStack
biomod.selected.models = bccvl.params$selected_models #'all' when all models have to be used to render projections or a subset vector of modeling.output models computed (eg, = grep('_RF', getModelsBuiltModels(myBiomodModelOut)))
# EMG If running one model at a time, this parameter becomes irrevelant
biomod.binary.meth = NULL #a vector of a subset of models evaluation method computed in model creation
biomod.filtered.meth = NULL #a vector of a subset of models evaluation method computed in model creation
biomod.compress = bccvl.params$compress # default 'gzip'; compression format of objects stored on your hard drive. May be one of `xz', `gzip' or NULL
biomod.build.clamping.mask = FALSE #if TRUE, a clamping mask will be saved on hard drive
opt.biomod.silent = FALSE #logical, if TRUE, console outputs are turned off
opt.biomod.do.stack = TRUE #logical, if TRUE, attempt to save all projections in a unique object i.e RasterStack
opt.biomod.keep.in.memory = TRUE #logical, if FALSE only the link pointing to a hard drive copy of projections are stored in output object
opt.biomod.output.format = NULL #'.Rdata', '.grd' or '.img'; if NULL, and new.env is not a Raster class, output is .RData defining projections saving format (on hard drive)


# model accuracy statistics
# these are available from dismo::evaluate.R NOT originally implemented in biomod2::Evaluate.models.R
dismo.eval.method = c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR")

# model accuracy statistics - combine stats from dismo and biomod2 for consistent output
model.accuracy = c(dismo.eval.method, biomod.models.eval.meth)
# TODO: these functions are used to evaluate the model ... configurable?

# read current climate data
current.climate.scenario = bccvl.enviro.stack(enviro.data.current, enviro.data.type, enviro.data.layer, resamplingflag=enviro.data.resampling)

###read in the necessary observation, background and environmental data
occur = bccvl.species.read(occur.data, month.filter) #read in the observation data lon/lat
absen = bccvl.species.read(absen.data, month.filter) #read in the observation data lon/lat

# geographically constrained modelling
if (!is.null(enviro.data.constraints) || enviro.data.generateCHall) {
  constrainedResults = bccvl.sdm.geoconstrained(current.climate.scenario, occur, absen, enviro.data.constraints, enviro.data.generateCHall);

  # Save a copy of the climate dataset
  current.climate.scenario.orig <- current.climate.scenario
  current.climate.scenario <- constrainedResults$raster
  occur <- constrainedResults$occur
  absen <- constrainedResults$absen
}

# Get the number of background points
pa_number_point = bccvl.params$pa_maxent_background_points

###run the models and store models
############### BIOMOD2 Models ###############
# 1. Format the data
# 2. Define the model options
# 3. Compute the model
# NOTE: Model evaluation is included as part of model creation

# BIOMOD_Modeling(data, models = c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT'), models.options = NULL,
#	NbRunEval=1, DataSplit=100, Yweights=NULL, Prevalence=NULL, VarImport=0, models.eval.meth = c('KAPPA','TSS','ROC'),
#	SaveObj = TRUE, rescal.all.models = TRUE, do.full.models = TRUE, modeling.id = as.character(format(Sys.time(), '%s')),
#	...)
#
# data	BIOMOD.formated.data object returned by BIOMOD_FormatingData
# models vector of models names choosen among 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF' and 'MAXENT'
# models.options BIOMOD.models.options object returned by BIOMOD_ModelingOptions
# NbRunEval	Number of Evaluation run
# DataSplit	% of data used to calibrate the models, the remaining part will be used for testing
# Yweights response points weights
# Prevalence either NULL (default) or a 0-1 numeric used to build 'weighted response weights'
# VarImport	Number of permutation to estimate variable importance
# models.eval.meth vector of names of evaluation metric among 'KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI' and 'ETS'
# SaveObj keep all results and outputs on hard drive or not (NOTE: strongly recommended)
# rescal.all.models	if true, all model prediction will be scaled with a binomial GLM
# do.full.models if true, models calibrated and evaluated with the whole dataset are done
# modeling.id character, the ID (=name) of modeling procedure. A random number by default.
# ... further arguments :
# DataSplitTable : a matrix, data.frame or a 3D array filled with TRUE/FALSE to specify which part of data must be used for models calibration (TRUE) and for models validation (FALSE). Each column correspund to a 'RUN'. If filled, args NbRunEval, DataSplit and do.full.models will be ignored.

###############
#
# Maxent
#
###############

# 1. Format the data as required by the biomod package
model.data = bccvl.biomod2.formatData(true.absen     = absen,
                                      pseudo.absen.points    = pa_number_point,
                                      pseudo.absen.strategy  = bccvl.params$pa_strategy,
                                      pseudo.absen.disk.min  = bccvl.params$pa_disk_min,
                                      pseudo.absen.disk.max  = bccvl.params$pa_disk_max,
                                      pseudo.absen.sre.quant = bccvl.params$pa_sre_quant,
                                      climate.data           = current.climate.scenario,
                                      occur                  = occur,
                                      species.name           = biomod.species.name,
                                      generate.background.data = TRUE,                # Generate background data as pseudo absence data
                                      species_algo_str       = species_algo_str)

# 2. Define the model options
model.options <- BIOMOD_ModelingOptions(MAXENT.Phillips = model.options.maxent)
# 3. Compute the model
model.sdm <-
    BIOMOD_Modeling(data = model.data,
                    models=c('MAXENT.Phillips'),
                    models.options=model.options,
                    NbRunEval=biomod.NbRunEval,
                    DataSplit=biomod.DataSplit,
                    Yweights=biomod.Yweights,
                    Prevalence=biomod.Prevalence,
                    VarImport=biomod.VarImport,
                    models.eval.meth=biomod.models.eval.meth,
                    SaveObj=TRUE,
                    rescal.all.models = biomod.rescal.all.models,
                    do.full.models = biomod.do.full.models,
                    modeling.id = biomod.modeling.id
                    )

# save the VIP plot
x.data <- attr(model.data,"data.env.var")
y.data <- attr(model.data,"data.species")
data1 = data.frame(y.data,x.data)
bccvl.VIPplot(method="maxent", data1=data1, pdf=TRUE,
              filename=paste('vip_plot', species_algo_str, sep="_"),
              this.dir=paste(biomod.species.name, "/models/bccvl", sep=""))

# model output saved as part of BIOMOD_Modeling() # EMG not sure how to retrieve
#save out the model object
bccvl.save(model.sdm, name="model.object.RData")

# Do projection over current climate scenario without constraint only if all env data layers are continuous.
if (enviro.data.genUnconstraintMap &&
   all(enviro.data.type == 'continuous') && 
   (!is.null(enviro.data.constraints) || enviro.data.generateCHall)) {
    model.proj <-
        BIOMOD_Projection(modeling.output     = model.sdm,
                          new.env             = current.climate.scenario.orig,
                          proj.name           = projection.name,
                          xy.new.env          = biomod.xy.new.env,
                          selected.models     = biomod.selected.models,
                          binary.meth         = biomod.binary.meth,
                          filtered.meth       = biomod.filtered.meth,
                          #compress            = biomod.compress,
                          build.clamping.mask = biomod.build.clamping.mask,
                          silent              = opt.biomod.silent,
                          do.stack            = opt.biomod.do.stack,
                          keep.in.memory      = opt.biomod.keep.in.memory,
                          output.format       = opt.biomod.output.format,
                          on_0_1000           = FALSE)

    # remove the current.climate.scenario to release disk space
    bccvl.remove.rasterObject(current.climate.scenario.orig)

    # convert projection output from grd to gtiff
    bccvl.grdtogtiff(file.path(getwd(),
                               biomod.species.name,
                               paste("proj", projection.name, sep="_")),
                     algorithm=ifelse(is.null(bccvl.params$subset), "maxent", sprintf("maxent_%s", bccvl.params$subset)),
                     filename_ext="unconstrained")

    # save the projection
    bccvl.saveProjection(model.proj, species_algo_str, filename_ext="unconstrained")
}

# predict for current climate scenario
model.proj <-
    BIOMOD_Projection(modeling.output=model.sdm,
                      new.env=current.climate.scenario,
                      proj.name  = projection.name,  #basename(enviro.data.current), {{ species }}
                      xy.new.env = biomod.xy.new.env,
                      selected.models = biomod.selected.models,
                      binary.meth = biomod.binary.meth,
                      filtered.meth = biomod.filtered.meth,
                      #compress = biomod.compress,
                      build.clamping.mask = biomod.build.clamping.mask,
                      silent = opt.biomod.silent,
                      do.stack = opt.biomod.do.stack,
                      keep.in.memory = opt.biomod.keep.in.memory,
                      output.format = opt.biomod.output.format,
                      on_0_1000 = FALSE)

# remove the current.climate.scenario to release disk space
bccvl.remove.rasterObject(current.climate.scenario)

# convert projection output from grd to gtiff
bccvl.grdtogtiff(file.path(getwd(),
                           biomod.species.name,
                           paste("proj", projection.name, sep="_")),
                 algorithm=ifelse(is.null(bccvl.params$subset), "maxent", sprintf("maxent_%s", bccvl.params$subset)))


# output is saved as part of the projection, format specified in arg 'opt.biomod.output.format'
loaded.model = BIOMOD_LoadModels(model.sdm, models="MAXENT.Phillips")
bccvl.saveBIOMODModelEvaluation(loaded.model, model.sdm, species_algo_str) 	# save output

# save the projection
bccvl.saveProjection(model.proj, species_algo_str)
