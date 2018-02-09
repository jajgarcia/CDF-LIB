#!/usr/bin/python
from math import *
import random


# ------------------------------------------------------------------------------

# read file with single column of data
#
# filename - name of input file
#
# return - array of data values
def read_dat1(filename):
  out=[]
  fin=open(filename,'r')
  for line in fin:
    if line[0] == "#": continue
    if line.strip() == "": continue
    tmp=line.split()
    out.append(float(tmp[0]))
  fin.close()
  return out

# ------------------------------------------------------------------------------

# write file with single column of data
#
# filename - name of input file
# dat - array with data points
#
# return - nothing
def write_dat1(filename,dat):
  nx=len(dat)
  fout=open(filename,'w')
  fout.write("#%10i\n"%nx)
  for x in dat:
    fout.write("%16.6E\n"%x)
  fout.close()

# ------------------------------------------------------------------------------

# read file with two columns of data
#
# filename - name of input file
#
# return - tuple of data value arrays
def read_dat2(filename):
  out1=[]
  out2=[]
  fin=open(filename,'r')
  for line in fin:
    if line[0] == "#": continue
    if line.strip() == "": continue
    tmp=line.split()
    out1.append(float(tmp[0]))
    out2.append(float(tmp[1]))
  fin.close()
  return (out1,out2)

# ------------------------------------------------------------------------------

# write file with double columns of data
#
# filename - name of input file
# dat1 - array with data to go in 1st column
# dat2 - array with data to go in 2nd column
#
# return - nothing
def write_dat2(filename,dat1,dat2):
  nx=len(dat1)
  if nx != len(dat2):
    print "data arrays must be of same length"
    raise
  fout=open(filename,'w')
  fout.write("#%10i\n"%nx)
  for ix in range(nx):
    fout.write("%16.6E%16.6E\n"%(dat1[ix],dat2[ix]))
  fout.close()

# ------------------------------------------------------------------------------

# create mesh with fixed spacing
#
# x0 - first point
# x1 - last point
# dx - mesh spacing
#
# return - array with mesh points
def make_mesh(x0,x1,dx):
  nx=int((x1-x0)/dx)+1
  mesh=[x0+dx*ix for ix in range(nx)]
  return mesh

# ------------------------------------------------------------------------------

# perform linear interpolation (or extrapolation)
#
# x0 - where to interpolate
# x1 - 1st x point
# y1 - 1st y point
# x2 - 2nd x point
# y2 - 2nd y point
#
# return - interpolated y point
def interpolate(x0,x1,y1,x2,y2):
  if (x1 == x2):
    print "cannot interpolate given indentical x1 and x2"
    raise
  slope=(y2-y1)/(x2-x1)
  return y1+slope*(x0-x1)

# ------------------------------------------------------------------------------

# create data set with constant profile; can optionally add result
# to existing profile (with base parameter)
# profile:  const
#
# const - profile value
# mesh - input mesh
# base - add distribution to this data (default: None)
#
# return - array with output profile
def pdf_constant(const,mesh,base=None):
  add=[0. for ix in range(len(mesh))]
  if base != None: add=base
  if (len(add) != len(mesh)):
    print "base distribution has wrong size"
    raise

  out=[a+const for a in add]
  return out

# ------------------------------------------------------------------------------

# create data set with Gaussian profile; can optionally add result
# to existing profile (with base parameter)
# Gaussian:  ac*exp(-(x-x0)**2/(2*ww**2))
#
# x0 - center of distribution
# ac - maximum value (at x=x0)
# ww - width
# mesh - input mesh
# base - add distribution to this data (default: None)
#
# return - array with output profile
def pdf_gaussian(x0,ww,ac,mesh,base=None):
  add=[0. for ix in range(len(mesh))]
  if base != None: add=base
  if (len(add) != len(mesh)):
    print "base distribution has wrong size"
    raise

  const=0.5/ww**2
  out=[add[ix]+ac*exp(-const*(mesh[ix]-x0)**2) for ix in range(len(mesh))]
  return out

# ------------------------------------------------------------------------------

# calculate Cumulative Distribution Function (CDF) of the given Probability
# Distribution Function (PDF)
# Note: assumes mesh spacing is fixed
#
# mesh - array with mesh points
# pdf - probability distribution function
#
# return - array with CDF
def calc_cdf(mesh,pdf):
  out=[0.5*(pdf[0]+pdf[1])]
  for ix in range(1,len(mesh)):
    out.append(out[ix-1]+0.5*(pdf[ix-1]+pdf[ix]))
  out=[y/out[-1] for y in out]
  return out

# ------------------------------------------------------------------------------

# lookup mesh value given probability
#
# prob - probability value: 0 - 1
# mesh - array with mesh points
# cdf - array with CDF
#
# return - data value associated with given probability
def cdf_lookup(prob,mesh,cdf):
  if prob < 0 or prob > 1:
    print "given probability, "+prob+", out of range"
    raise

  # due to CDF integration accuracy, given probability can sometimes
  # be larger than the last CDF point (but still less than 1); in this
  # case, return the largest mesh point
  if prob > cdf[-1]: return mesh[-1]

  ix=0
  while cdf[ix] < prob: ix+=1
  out=interpolate(prob,cdf[ix-1],mesh[ix-1],cdf[ix],mesh[ix])
#  slope=(mesh[ix]-mesh[ix-1])/(cdf[ix]-cdf[ix-1])
#  out=mesh[ix-1]+slope*(prob-cdf[ix-1])
  if out < mesh[0] or out > mesh[-1]: out=None
  return out

# ------------------------------------------------------------------------------

# generate list of random events using given CDF
# 
# mesh - array with mesh points
# cdf - array with CDF
# npoints - number of events to generate
#
# return - list of events
def gen_events(mesh,cdf,npoints):
  out=[]
  n=0
  while 1:
    if n > npoints: break
    r=random.random()
    v=cdf_lookup(r,mesh,cdf)
    if v is None: continue      # lookup out-of-range
    n+=1
    out.append(v)
  return out

# ------------------------------------------------------------------------------

# bin event data using given mesh
# note: assumes mesh has constant spacing
#
# mesh - array with mesh points
# events - list of events
#
# return - tuple(bin mesh, binned data)
def bin_events(mesh,events):
  events.sort()

  # construct bin mesh
  bmesh=[0.5*(mesh[ix]+mesh[ix+1]) for ix in range(len(mesh)-1)]

  out=[0. for ib in range(len(bmesh))]
  ix=0
  for e in events:
    while e > mesh[ix]: ix+=1
    out[ix-1]+=1
  return (bmesh,out)

# ------------------------------------------------------------------------------

# get average from event list
#
# events - list of events
#
# return - average value
def avg_events(events):
  return sum(events)/float(len(events))


# ------------------------------------------------------------------------------

# get average from PDF (or binned event list)
#
# mesh - array of mesh points
# pdf - list of events
#
# return - average value
def avg_pdf(mesh,pdf):
  tmp=[mesh[ix]*pdf[ix] for ix in range(len(mesh))]
  return sum(tmp)/float(sum(pdf))

# ------------------------------------------------------------------------------

# compute 3-point derivative of given function (assumes mesh with fixed
# spacing)
#
# mesh - array of mesh points
# f - array of function points
#
# return - array of derivative points
def der1_3pt(mesh,f):
  out=[]
  dx=mesh[1]-mesh[0]
  factor=0.5/dx

  # first point uses end-point formula
  out.append(factor*(-3.*f[0]+4.*f[1]-f[2]))

  # middle points (not 1st or last)
  for ix in range(1,len(mesh)-1):
    out.append(factor*(f[ix+1]-f[ix-1]))

  # last point uses end-point formula
  out.append(factor*(3.*f[-1]-4.*f[-2]+f[-3]))

  return out

# ------------------------------------------------------------------------------

# shift function by given delta-x
#
# mesh - array of mesh points
# pdf - array with PDF
# dx - amount of shift
#
# return - array of shifted points
def shift_dat(mesh,pdf,dx):
  out=[]
  for ix in range(len(mesh)):
    x=mesh[ix]
    xp=x-dx

    if xp < mesh[0]:           # constant background to left
      out.append(pdf[0])
    elif xp > mesh[-1]:        # constant background to right
      out.append(pdf[-1])
    else:                      # normal interpolation
      ixp=0
      while mesh[ixp] < xp: ixp+=1
      out.append(interpolate(xp,mesh[ixp-1],pdf[ixp-1],mesh[ixp],pdf[ixp]))

  return out

# ------------------------------------------------------------------------------

# regrid data set
#
# mesh - original mesh
# dat - original data
# pmesh - new mesh
#
# return - array of data on new mesh
def regrid_dat(mesh,dat,pmesh):
  out=[]
  for ix in range(len(pmesh)):
    x=pmesh[ix]

    if x < mesh[0]:            # interpolate to the left
      out.append(interpolate(x,mesh[0],dat[0],mesh[1],dat[1]))
    elif x > mesh[-1]:         # interpolate to the right
      out.append(interpolate(x,mesh[-2],dat[-2],mesh[-1],dat[-1]))
    else:                      # normal interpolation
      ixp=0
      while mesh[ixp] < x: ixp+=1
      out.append(interpolate(x,mesh[ixp-1],dat[ixp-1],mesh[ixp],dat[ixp]))

  return out

# ------------------------------------------------------------------------------

# calculate delta-x functional for fitting profile to data
# also return the profile scale parameter and chi^2
#
# mesh - array with mesh points
# g - data points
# f - profile points
# df - 1st derivative of profile
#
# return - tuple: (value of functional, scale factor, chi^2)
def func_dx(mesh,g,f,df):

  sum_gf=0.
  sum_fdf=0.
  sum_gdf=0.
  sum_ff=0.
  for ix in range(len(mesh)):
    sum_gf+=g[ix]*f[ix]
    sum_fdf+=f[ix]*df[ix]
    sum_gdf+=g[ix]*df[ix]
    sum_ff+=f[ix]*f[ix]

  z=sum_gf*sum_fdf-sum_gdf*sum_ff
  scla=sum_gf/sum_ff
  return (z,scla)

# ------------------------------------------------------------------------------

# use secant method to fit a profile to data set
#
# mesh - array with mesh points
# pdf - profile
# dat - data set
# dx1 - 1st initial condition of profile shift
# dx2 - 2nd initial condition of profile shift
#
# return - tuple: (fitted shift, fitted scale factor, fitted pdf)
def fit_pdf(mesh,pdf,dat,dx1,dx2):

  # set tolerance 
  # method stops once relative change in functional is less than
  # this amount
  tol=0.01
  maxit=20    # maximum number of iterations

  # calculate fitting functional for 1st initial condition
  pdfp=shift_dat(mesh,pdf,dx1)
  dpdf=der1_3pt(mesh,pdfp)
  (z1,scl1)=func_dx(mesh,dat,pdfp,dpdf)

  # calculate fitting functional for 2nd initial condition
  pdfp=shift_dat(mesh,pdf,dx2)
  dpdf=der1_3pt(mesh,pdfp)
  (z2,scl2)=func_dx(mesh,dat,pdfp,dpdf)

  # start secant loop
  it=0
  while 1:

    # not converging?
    it+=1
    if (it > maxit):
      print "secant method does not converge"
      raise

    # calculate relative change in functional to check if finished
    dz=2.*abs(z2-z1)/abs(z2+z1)
    if dz <= tol: break

    # get next dx    
    dxp=dx2-z2*(dx2-dx1)/(z2-z1)

    # shift parameters
    dx1=dx2
    z1=z2
    scl1=scl2 

    # calculate next functional
    dx2=dxp
    pdfp=shift_dat(mesh,pdf,dx2)
    dpdf=der1_3pt(mesh,pdfp)
    (z2,scl2)=func_dx(mesh,dat,pdfp,dpdf)

  print ""
  print "----------------------------------"
  print "secant method finished: "
  print "  num iterations: %5i"%it
  print "  functional: %16.6E"%z2
  print "  relative change: %12.2E"%dz
  print "----------------------------------"
  print ""

  # apply scale to fitted profile
  pdfp=[scl2*y for y in pdfp]

  return (dx2,scl2,pdfp)

# ------------------------------------------------------------------------------

# calculate chi^2
#
# mesh - array with mesh points
# pdf - profile
# dat - data set
#
# return - reduced chi2
def calc_chi2(mesh,pdf,dat):
  sum_num=0.
  sum_obs=0.
  sum_obs2=0.
  for ix in range(len(mesh)):
    sum_num+=(dat[ix]-pdf[ix])**2
    sum_obs+=dat[ix]
    sum_obs2+=dat[ix]**2
  nfac=1./float(len(mesh))
  var=nfac*sum_obs2-nfac**2*sum_obs**2
  chi2=sum_num/var
  return chi2

# ------------------------------------------------------------------------------

# calculate R^2
#
# mesh - array with mesh points
# pdf - profile
# dat - data set
#
# return - reduced chi2
def calc_r2(mesh,pdf,dat):
  avgdat=sum(dat)/float(len(mesh))
  sum_err=0.
  sum_tot=0.
  for ix in range(len(mesh)):
    sum_err+=(dat[ix]-pdf[ix])**2
    sum_tot+=(dat[ix]-avgdat)**2
  return 1.-sum_err/sum_tot







