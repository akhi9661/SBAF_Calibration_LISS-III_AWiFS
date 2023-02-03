'''

This code module applies a "simple Spectral Band Adjustment Factor (SBAF)" on LISS III and AWiFS to calibrate them with the help
of a reference image like Landsat 8 and Sentinel 2.

Inputs:

inpf_liss = path to folder LISS III or AWiFS bands
inpf_ref = path to folder containing reference image bands

Output:
Final output folder: Georeferenced


Note: Please make sure that both input and reference image have the same projection system. 

'''

inpf_liss = r'D:\Aerosol Modelling\Aerosol\Output\LISS III\226922421'
inpf_ref = r'D:\Aerosol Modelling\Aerosol\Output\S2A1C\T43RGM_Delhi_2021\T43RGM_2021-11-25_Delhi\Test'

''' ------------------------------------------------------------------------------------------------------------------ '''

def meta(inpf, keyword):
    
    file = open(glob.glob(os.path.join(inpf, '*_META.txt'))[0]).readlines()
    meta = ""
    
    for lines in file:
        if keyword in lines:
            meta = float(lines.split("=")[-1].strip())
    return meta

def toa_reflect(inpf, inp_name, opf, band_no):
    
    esol = {'B2': 1849.5, 'B3': 1553.0, 'B4': 1092.0, 'B5': 239.52}
    esol_band = list(esol.values())[band_no-2]
    
    lmax, lmin, sun_elev = meta(inpf, f'B{band_no}_Lmax'), meta(inpf, f'B{band_no}_Lmin'), meta(inpf, 'SunElevationAtCenter')
    
    with rasterio.open(os.path.join(inpf, inp_name)) as r:
        param = r.profile
        toa_raw = r.read(1).astype('float32')
    
    param.update(dtype = 'float32')
    toa_raw[toa_raw == 0] = np.nan
    toa_rad = lmin + ((lmax - lmin)/1024)*toa_raw
    reflectance = (np.pi * 1 * toa_rad) / (esol_band * sin(radians(sun_elev)))
    
    reflectance[reflectance>1] = np.nan
    reflectance[reflectance<0] = np.nan
    if (np.nanmax(reflectance) <= np.nanpercentile(reflectance, 99.99)):
        reflectance = reflectance
    else:
        reflectance[reflectance>=np.nanpercentile(reflectance, 99.99)] = np.nanpercentile(reflectance, 99.999)
    
    op_name = os.path.basename(inp_name).split('.')[0] + '_ref.TIF'
    with rasterio.open(os.path.join(opf, op_name), 'w', **param) as r:
        r.write(reflectance, 1)
    
    return 'Done'

def do_ref(inpf, opf):
    
    print('Radiance to reflectance conversion: LISS III.')
    
    original = os.listdir(inpf)
    gtif = list(filter(lambda x: x.endswith(("tif", "TIF", "img")), original))
    for gi in gtif:
        band_no = int(''.join(list(filter(str.isdigit, gi))))
        toa_reflect(inpf, gi, opf, band_no)
    
    return None
    
def create_multiband_image(inpf_liss, inpf_ref, files_liss, files_ref):
    
    ''' LISS III Composite '''
    
    print('Stacking: LISS III.')
    with rasterio.open(os.path.join(inpf_liss, files_liss[0])) as src:
        profile = src.profile
        multi_band_liss = np.zeros((profile['height'], profile['width'], len(files_liss)))

    for i, filename in enumerate(files_liss):
        with rasterio.open(os.path.join(inpf_liss, filename)) as src:
            multi_band_liss[:,:,i] = src.read(1)
            
    profile.update(count = multi_band_liss.shape[2])
    profile.update(dtype = 'float32')
    op_liss = os.path.join(inpf_liss, 'composite_liss.TIF')
    with rasterio.open(op_liss, 'w', **profile) as dst:
        dst.write(np.rollaxis(multi_band_liss.astype(profile['dtype']), axis=2))
    dst.close()
        
    ''' Reference Image Composite '''
    
    print('Stacking: Reference image.')
    with rasterio.open(os.path.join(inpf_ref, files_ref[0])) as src:
        profile = src.profile
        multi_band_ref = np.zeros((profile['height'], profile['width'], len(files_ref)))

    for i, filename in enumerate(files_ref):
        with rasterio.open(os.path.join(inpf_ref, filename)) as src:
            multi_band_ref[:,:,i] = src.read(1).astype('float32')*0.0001
            
    profile.update(count = multi_band_ref.shape[2])
    profile.update(dtype = 'float32')
    op_ref = os.path.join(inpf_ref, 'composite_ref.TIF')
    with rasterio.open(op_ref, 'w', **profile) as dst:
        dst.write(np.rollaxis(multi_band_ref, axis=2))
    dst.close()
    
    return (op_liss, op_ref)

def do_multiband(inpf_liss, inpf_ref):
    
    original = os.listdir(inpf_liss)
    gtif_liss = list(filter(lambda x: x.endswith(("tif", "TIF", "img")), original))
    
    original = os.listdir(inpf_ref)
    gtif_ref = list(filter(lambda x: x.endswith(("tif", "TIF", "img")), original))
    
    op_liss, op_ref = create_multiband_image(inpf_liss, inpf_ref, gtif_liss, gtif_ref)
    return (op_liss, op_ref)
    
def resample_image(file_liss, file_ref):
    
    print('Resampling: LISS III to reference image')
    
    reference = gdal.Open(file_ref, 0)
    referenceTrans = reference.GetGeoTransform()
    x_res = referenceTrans[1]
    y_res = -referenceTrans[5]

    opf_resample = os.path.join(os.path.dirname(file_liss), os.path.basename(file_liss).split('.')[0] + '_resample.TIF')

    kwargs = {"format": "GTiff", "xRes": x_res, "yRes": y_res}
    ds = gdal.Warp(opf_resample, file_liss, **kwargs)
    ds = None
    
    return opf_resample

def clip_image(opf_resample, file_ref):
    
    print('Clipping: LISS III to reference image')
    
    maskDs = gdal.Open(file_ref, 0)
    projection = maskDs.GetProjectionRef()
    geoTransform = maskDs.GetGeoTransform()
    minx = geoTransform[0]
    maxy = geoTransform[3]
    maxx = minx + geoTransform[1] * maskDs.RasterXSize
    miny = maxy + geoTransform[5] * maskDs.RasterYSize


    data = gdal.Open(opf_resample, 0)
    opf_clip = os.path.join(os.path.dirname(opf_resample), os.path.basename(opf_resample).split('.')[0] + '_clip.TIF')
    ds = gdal.Translate(opf_clip, data, format = 'GTiff',projWin = [minx,maxy,maxx,miny], outputSRS = projection) 
    ds = None
    
    return opf_clip
    

def calc_sbaf(file_liss, file_ref):
    
    opf_resample = resample_image(file_liss, file_ref)
    opf_clip = clip_image(opf_resample, file_ref)
    
    print('Calculating SBAF.')

    with rasterio.open(opf_clip) as r:
        liss = r.read().astype('float32')
        param = r.profile
        print('LISS III shape:', liss.shape)
    with rasterio.open(file_ref) as r:
        ref = r.read().astype('float32')
        print('Reference shape:', ref.shape)
        
    cal_path = os.path.join(os.path.dirname(file_liss), 'Calibrated')
    if os.path.exists(cal_path):
        shutil.rmtree(cal_path)
    os.makedirs(cal_path)
        
    num_bands = liss.shape[0]
    cal_liss = np.zeros((liss.shape[1], liss.shape[2]), dtype = 'float32')
    sbaf = []
    for i in range(num_bands):
        band_liss = liss[i,:,:]
        band_reference = ref[i,:,:]
        sbaf.append(np.nanmean(band_reference)/np.nanmean(band_liss))
        sbaf_temp = sbaf[i]
        
        print(f'Calibrating band {i+2}...')
        print('Factor:', sbaf_temp)
        cal_liss = band_liss * sbaf_temp
        opf_cal = os.path.join(cal_path, f'Band_{i+2}_cal.TIF')
        with rasterio.open(opf_cal, 'w', **param) as r:
            r.write(cal_liss, 1)

    
    os.remove(opf_clip)
    os.remove(opf_resample)
    os.remove(file_liss)
    os.remove(file_ref)
    
    print("'Done'")
    return (sbaf, cal_liss, band_liss, band_reference)


opf = os.path.join(inpf_liss, 'Reflectance')
if os.path.exists(opf):
    shutil.rmtree(opf)
os.makedirs(opf)

do_ref(inpf_liss, opf)
op_liss, op_ref = do_multiband(opf, inpf_ref)
sbaf, cal_liss, ref_liss, ref_band = calc_sbaf(op_liss, op_ref)