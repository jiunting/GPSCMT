#!/usr/bin/python
import math
import numpy




def kagan_angle(str1,dip1,rak1,str2,dip2,rak2):
    str1 = math.radians(str1)
    dip1 = math.radians(dip1)
    rak1 = math.radians(rak1)
    str2 = math.radians(str2)
    dip2 = math.radians(dip2)
    rak2 = math.radians(rak2)
    n1 = numpy.array([-math.sin(dip1)*math.sin(str1), -math.sin(dip1)*math.cos(str1), math.cos(dip1)])
    n2 = numpy.array([-math.sin(dip2)*math.sin(str2), -math.sin(dip2)*math.cos(str2), math.cos(dip2)])
    d1 = numpy.array([math.cos(rak1)*math.cos(str1)+math.sin(rak1)*math.cos(dip1)*math.sin(str1), -math.cos(rak1)*math.sin(str1)+math.sin(rak1)*math.cos(dip1)*math.cos(str1), math.sin(rak1)*math.sin(dip1)])
    d2 = numpy.array([math.cos(rak2)*math.cos(str2)+math.sin(rak2)*math.cos(dip2)*math.sin(str2), -math.cos(rak2)*math.sin(str2)+math.sin(rak2)*math.cos(dip2)*math.cos(str2), math.sin(rak2)*math.sin(dip2)])
    MT1 = numpy.array([[2*n1[0]*d1[0], n1[0]*d1[1]+n1[1]*d1[0], n1[0]*d1[2]+n1[2]*d1[0]],[n1[0]*d1[1]+n1[1]*d1[0],2*n1[1]*d1[1],n1[1]*d1[2]+n1[2]*d1[1]],[n1[0]*d1[2]+n1[2]*d1[0],n1[1]*d1[2]+n1[2]*d1[1],2*n1[2]*d1[2]]])
    MT2 = numpy.array([[2*n2[0]*d2[0], n2[0]*d2[1]+n2[1]*d2[0], n2[0]*d2[2]+n2[2]*d2[0]],[n2[0]*d2[1]+n2[1]*d2[0],2*n2[1]*d2[1],n2[1]*d2[2]+n2[2]*d2[1]],[n2[0]*d2[2]+n2[2]*d2[0],n2[1]*d2[2]+n2[2]*d2[1],2*n2[2]*d2[2]]])
    cosom = numpy.sum(numpy.multiply(MT1,MT2))/numpy.linalg.norm(MT1)/numpy.linalg.norm(MT2)
    omega = math.degrees(math.acos(cosom))
    return(omega)

