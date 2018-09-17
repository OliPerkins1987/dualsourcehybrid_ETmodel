# dualsourcehybrid_ETmodel
This is a PC raster implementation of the dual source hybrid (Penman Monteith) method for modelling ET (EG Guan and Wilson 2009). 

To install pcraster please go to http://pcraster.geo.uu.nl/

The model requires LULC and DEM rasters as inputs, and is forced with wide range of meteorological inputs. Comments in the model point towards relevant academic literature for sources. In order to function correctly all model inputs need to be stored in the same folder on your hard drive. 

The model is intended for use in montane environments, and places significant emphasis on the imapct of topography. A current limitation is that runoff is immediately removed from the model system, rather than flowing to adjacent raster cells. 

Please contact me with any questions, comments or suggestions at oliver.perkins@kcl.ac.uk




