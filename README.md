# Mounty: a PCraster dual source hybrid Penman-Monteith model of ET, designed for use in montane environments

This is a PCraster implementation of the dual source hybrid (Penman Monteith) method for modelling ET (EG Guan and Wilson 2009). It is a coupled model, that combines simple representations of the atmosphere, soil and vegetation. An overview of model function and conceptualisation is given in two simple flow diagrams. 

The model requires LULC and DEM rasters as inputs. It attempts to adopt a parsimonious approach to metereological forcing requirements, and uses the dew point method to avoid estimates of Relative Humidity (See Allen et al., 1998). It requires

Comments in the model point towards relevant academic literature for sources. 

The model is intended for use in montane environments, and places significant emphasis on the imapct of topography. A current limitation is that runoff is immediately removed from the model system, rather than flowing to adjacent raster cells. 

To install pcraster please go to http://pcraster.geo.uu.nl/ 

In order to function correctly all model inputs need to be stored in the same folder on your hard drive. You can download example data to force the model, to demonstrate model function and correct formatting for PCRaster from https://www.dropbox.com/sh/s2yqshcqgokssb0/AAAVxJ8xffH80kuZ9FwEyZOCa?dl=0.


Please contact me with any questions, comments or suggestions at oliver.perkins@kcl.ac.uk




