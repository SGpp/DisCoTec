import Tkinter
import math

class SimplePlot(Tkinter.Canvas):

    def plot(self, x, y,x0,y0,xsize,ysize,maxval,minval):
        for i in range(len(x)-1):
            if (maxval!=0 and y[i]>=minval):
                self.create_line(x0+xsize*x[i]/max(x), y0+ysize-ysize*(y[i]-minval)/(maxval-minval), 
                                 x0+xsize*x[i+1]/max(x), y0+ysize-ysize*(y[i+1]-minval)/(maxval-minval), fill="red", width=2)
            else:
                self.create_line(x0+xsize*x[i]/max(x), y0+ysize, 
                                 x0+xsize*x[i+1]/max(x), y0+ysize, fill="red", width=2)

    def create_axis(self,x0,y0,xsize,ysize,x,y,quant,profgrad,log):
        num_tags_x=5
        maxtagx=1.0
        for i in range(num_tags_x+1):
            posx=x0+xsize*i/num_tags_x*maxtagx
            posy=y0+ysize+5
            self.create_text((posx,posy),anchor='n',text=str(maxtagx*i/num_tags_x))
            self.create_line(posx,y0+ysize,posx,y0)            
        num_tags_y=5
        fac=1.0
        if not log:
            if max(y)>0.:
                while max(y)*fac>20:
                    fac=fac/10
                while max(y)*fac<2:
                    fac=fac*10
            else:
                while max(y)*fac<-20:
                    fac=fac/10
                while max(y)*fac>-2:
                    fac=fac*10                
        maxtagy=round(fac*max(y))
        if log:
            minval=max(min(y),-2)
        else:
            minval=min(0,min(y))
        mintagy=round(fac*minval)
        minval=mintagy/fac
        maxval=max(max(y),maxtagy/fac)
        difftags=(maxtagy/fac-mintagy/fac)/(num_tags_y-1.0)
        for i in range(num_tags_y):
            posx=x0-5
            if max(y)!=0:
                posy=ysize+y0-float(i)*ysize*difftags/(maxval-minval)
            else:
                posy=ysize+y0
            if log:
                self.create_text((posx,posy),anchor='e',text=str('%1.2f' % math.pow(10.,mintagy/fac+difftags*i)))
            else:
                self.create_text((posx,posy),anchor='e',text=str('%1.1f' % (mintagy/fac+difftags*i)))
            self.create_line(x0,posy,x0+xsize,posy)
        if quant.get()==0:
            title='Temperature'
        elif quant.get()==1:
            title='Density'
        else:
            title='Heat source profile'
        if quant.get()!=2:
            if profgrad==0:
                title+=' profile'
            elif profgrad==1:
                title+=' gradient'
        self.create_text((x0+xsize/2.,y0/2),text=title)
        self.create_text((x0+xsize/2.,y0+ysize+25),text='r/a')
        #self.create_text((x0/2,y0+ysize/2)

        return maxval,minval

