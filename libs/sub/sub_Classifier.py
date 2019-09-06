from numpy import array,column_stack
from numpy.random import randint

def get_random_time_segment(segment_frames,total_frames=12000):   
    '''
    Gets a random time segment of duration segment_frames in a file
    with number of frames: total_frames
    '''
    
    segment_start = randint(0, high = total_frames-
                                   segment_frames)
    segment_end = segment_start + segment_frames
    
    return (segment_start, segment_end)


def selROI(movie,x,y,t = 0,tw = None,xlim = (0,512),ylim= (0,256),size = 1):
    '''Select ROI from a movie with shape (time,height,width)
    '''
    if x<xlim[0]+size:
        x0,xf = (xlim[0],xlim[0]+2*size+1)
    elif x>xlim[1]-size-1:
        x0,xf = (xlim[1]-2*size-1,xlim[1])
    else:
        x0,xf = (x-size,x+size+1)
    if y<ylim[0]+size:
        y0,yf = (ylim[0],ylim[0]+2*size+1)
    elif y>ylim[1]-size-1:
        y0,yf = (ylim[1]-2*size-1,ylim[1])
    else:
        y0,yf = (y-size,y+size+1)
    
    if tw is None:
        return(movie[t:,y0:yf,x0:xf])
    else:
        return(movie[t:(t+tw),y0:yf,x0:xf])       

def extract_data_Class20181121(roi):
    ''' Provide data format to Classifier 20181121_17.h5'''
    
    start,end = get_random_time_segment(1000,total_frames=roi.shape[0])
    xtt = roi[start:end,:].transpose()
    s0 = xtt.std()/xtt.mean()
    s1 = xtt.mean()
    xtt = xtt/s1
    s2 = xtt.std()
    s3 = ((xtt-1.0)**3).mean()/s2**3
    s4 = ((xtt-1.0)**4).mean()/s2**4
    s0t = array([s2,s3,s4,(s3**2+1/s4),0,0,0,0,0])
    
    return(column_stack((s0t,xtt)))
    
def formatClass20181121(Xt):
    ''' Provide data format prior to classify to Classifier 20181121_17.h5'''
    Xt = array(Xt)
    Xt = Xt.reshape((Xt.shape[0],Xt.shape[1],Xt.shape[2],1))

    return(Xt)
