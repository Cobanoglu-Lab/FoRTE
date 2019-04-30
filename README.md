# FoRTE

### FunctiOnal, physiologically Relevant Therapeutic Evaluation

The software is intended to count pixels in multiple z-stacks of images. It is the companion analysis software for the high-content assay described in: https://www.biorxiv.org/content/10.1101/312504v3


Usage:

```python process_images.py <FOLDER> <OUTPUT_FOLDER> <CHN1> <CHN2>...```

where CHN1, CHN2, etc are the names of the channels. 

Example: 

```python process_images.py DMSO output Nuc Apop EtHd```


Inside the folder, the images should be organized as: 

```xy<XY_COORD>z<Z_COORD>c<CHANNEL_ID>.tif```

Example:

```xy01z001c1.tif```


When running the usage example with the example image names, ```c1``` gets mapped to the name ```Nuc```, ```c2``` gets mapped to ```Apop```, and ```c3``` to ```EtHd```.  The result files will be formatted appropriately. 

