# -*- coding: utf-8 -*-
"""
Created on Sun Aug 29 10:46:37 2021

@author: ENGINEERING
"""
import math
import numpy as np
import matplotlib.pyplot as plt
input_angle=0
output_angle=0
def input_angles_from_cheb(numberOfAngles=3,initial_input_angle=180, final_input_angle=120):
    cheb_input_angles=[0]*numberOfAngles
    for j in range(1,numberOfAngles+1):
        first_section=(initial_input_angle+final_input_angle)*0.5
        cos_part=((2*j-1)/6)*180
        #print(cos_part)
        half_difference=(final_input_angle-initial_input_angle)*0.5
        second_section=half_difference*math.cos(cos_part)
        angles=first_section-second_section
        cheb_input_angles[j-1]=angles
    print("cheb input angles:")
    print(cheb_input_angles)
    return cheb_input_angles

def Output_angles_from_cheb(numberOfAngles=3,initial_input_angle=180, final_input_angle=120):
    cheb_output_angles=[0]*numberOfAngles
    chebs=input_angles_from_cheb(numberOfAngles,initial_input_angle, final_input_angle)
    for i in range(1,numberOfAngles+1):
        cheb_output_angles[i-1]=70-(18000/chebs[i-1])
    #print("cheb output angles:")
    #print(cheb_output_angles)
    return cheb_output_angles


def calculateKs(numberOfAngles=3,initial_input_angle=180, final_input_angle=120):
    LHS=np.array([[0.00, 0.00, 0], [0.00, 0.00, 0], [0.00, 0.00, 0]])
    RHS=np.array([0.00, 0.00, 0.])
    input_angles=input_angles_from_cheb(numberOfAngles,initial_input_angle, final_input_angle)
    output_angles=Output_angles_from_cheb(numberOfAngles,initial_input_angle, final_input_angle)
    for i in range(1,4):
        #print(math.cos(input_angles[i-1]))
        #print(math.cos(output_angles[i-1]))
        #print(math.cos(math.radians(154)))
        LHS[i-1,0]=math.cos(math.radians(output_angles[i-1]))
        LHS[i-1,1]=-(math.cos(math.radians(input_angles[i-1])))
        LHS[i-1,2]=1
    
    for j in range(1,4):
        RHS[j-1]=math.cos(math.radians(input_angles[j-1]-output_angles[j-1]))
    #print(RHS)
    array_of_ks = np.linalg.solve(LHS, RHS)
    print("Ks From Normal Array Manipulation")
    print(array_of_ks)
    return array_of_ks



def calculate_link_lengths(numberOfAngles=3,initial_input_angle=180, final_input_angle=120, least_squareMethod=False):
    if (least_squareMethod==False):
        all_Ks=calculateKs(numberOfAngles,initial_input_angle, final_input_angle)
    else:
        all_Ks=least_square_method(numberOfAngles)
    length_array=[0]*4
    k1=all_Ks[0]
    k2=all_Ks[1]
    k3=all_Ks[2]
    a=600#mm
    length_array[0]=a
    d=a*all_Ks[0]
    length_array[1]=d
    c=abs(d*all_Ks[1])
    length_array[2]=c
    b=math.sqrt((a*a)+(c*c)+(d*d)-(k3*2*a*c))
    length_array[3]=d
    print ("Link Lengths[a,b,c,d]:")
    print(length_array)
    return length_array
def transmission_angle(numberOfAngles=3,initial_input_angle=180, final_input_angle=120):
    length_array=calculate_link_lengths(numberOfAngles,initial_input_angle, final_input_angle)
    a=length_array[0]
    b=length_array[1]
    c=length_array[2]
    d=length_array[3]
    asq=a*a
    bsq=b*b
    csq=c*c
    dsq=d*d
    pos=0
    Intervals=5
    array_size=(180-120)/Intervals
    input_Angle=[0]*int(array_size)
    transmision_Angle=[0]*int(array_size)
    for i in range (120,180,Intervals):
        input_Angle[pos]=i
        numerator=bsq+csq-asq-dsq+(2*a*d*math.cos(math.radians(i)))
        transmision_angle=numerator/(2*b*c)
        transmision_Angle[pos]=transmision_angle
        pos=pos+1
    print("Calculated Transmision Angles:")
    print(transmision_Angle) 
    plt.plot(input_Angle, transmision_Angle)
    plt.show()

def least_square_method(numberOfAngles=5,initial_input_angle=180, final_input_angle=120):
    inputAngles=input_angles_from_cheb(numberOfAngles,initial_input_angle, final_input_angle)
    outputAngles=Output_angles_from_cheb(numberOfAngles,initial_input_angle, final_input_angle)
    LHS=np.array([[0.00, 0.00, 0], [0.00, 0.00, 0], [0.00, 0.00, 0]])
    RHS=np.array([0.00, 0.00, 0.])
    sum_of_cossqT4, sum_ofcosT2T4, sum_ofcosT2, sum_ofcossqT2, sum_ofcosT4=0,0,0,0,0
    N=numberOfAngles
    sum_ofCT4CT2_T4, sum_ofCT2CT2_T4, sum_ofCT2_T4=0,0,0
    for i in range(1,numberOfAngles+1):
        sum_ofcosT2=sum_ofcosT2 +math.cos(math.radians(inputAngles[i-1])) #checked for accuracy
        sum_ofcosT4=sum_ofcosT4 +math.cos(math.radians(outputAngles[i-1])) #checked for accuracy
        sum_ofcosT2T4=sum_ofcosT2T4 +math.cos(math.radians(inputAngles[i-1]))*math.cos(math.radians(outputAngles[i-1]))
        sum_ofcossqT2=sum_ofcossqT2 +math.cos(math.radians(inputAngles[i-1]))**2 
        sum_of_cossqT4=sum_of_cossqT4 +math.cos(math.radians(outputAngles[i-1]))**2
        sum_ofCT4CT2_T4=sum_ofCT4CT2_T4 + math.cos(math.radians(outputAngles[i-1]))*math.cos(math.radians(inputAngles[i-1]-outputAngles[i-1]))
        sum_ofCT2CT2_T4=sum_ofCT2CT2_T4 + math.cos(math.radians(inputAngles[i-1]))*math.cos(math.radians(inputAngles[i-1]-outputAngles[i-1]))
        sum_ofCT2_T4=sum_ofCT2_T4 + math.cos(math.radians(inputAngles[i-1]-outputAngles[i-1]))
    
    LHS[0,0]=sum_of_cossqT4
    LHS[0,1]=-(sum_ofcosT2T4)
    LHS[0,2]=sum_ofcosT4
    LHS[1,0]=sum_ofcosT2T4
    LHS[1,1]=-(sum_ofcossqT2)
    LHS[1,2]=sum_ofcosT2
    LHS[2,0]=sum_ofcosT4
    LHS[2,1]=-(sum_ofcosT2)
    LHS[2,2]=N
    RHS[0]=sum_ofCT4CT2_T4
    RHS[1]=sum_ofCT2CT2_T4
    RHS[2]=sum_ofCT2_T4
    array_of_Ks = np.linalg.solve(LHS, RHS)
    print("Ks From Least Square Method")
    print(array_of_Ks)
    return array_of_Ks
def structural_errors(numberOfAngles=3,initial_input_angle=180, final_input_angle=120, Range=5):
    array_size=abs(int((final_input_angle-initial_input_angle)/Range))
    expected_output_angles=Required_Output_Angles(Range=Range)
    calculated_output_angles=Generated_Output_Angles(Range=Range,numberOfAngles=numberOfAngles)
    errors=[0.00]*array_size
    index_counter=0
    for i in expected_output_angles:
        errors[index_counter] =calculated_output_angles[index_counter]-expected_output_angles[index_counter]
        index_counter=index_counter+1
    print ("Structural Errors")
    print(errors)
    return errors
def Required_Output_Angles(Range=5,initial_input_angle=180, final_input_angle=120):
    array_size=abs(int((final_input_angle-initial_input_angle)/Range))
    expected_output_angles=[0.00]*array_size
    index_counter=0
    input_angles=input_angle_generator(Range,initial_input_angle,final_input_angle)
    for i in input_angles:
        expected_output_angles[index_counter]=70-(18000/i)
        index_counter=index_counter+1
    return expected_output_angles
def Generated_Output_Angles(numberOfAngles=3,Range=5,initial_input_angle=180, final_input_angle=120):
    if (numberOfAngles>3):
        Ks=least_square_method(numberOfAngles=numberOfAngles)
        
    else:
        Ks=calculateKs()
    k1=Ks[0]
    k2=Ks[1]
    k3=Ks[2]
    A=0
    B=0
    C=0
    array_size=abs(int((final_input_angle-initial_input_angle)/Range))
    calculated_output_angles=[0.00]*array_size
    index_counter=0
    expected_output_angles=Required_Output_Angles(Range=Range)
    for i in input_angle_generator(Range,initial_input_angle,final_input_angle):
        A=(1-k2)*math.cos(math.radians(i))-k1+k3#checked for accuracy
        B=-(2*math.sin(math.radians(i))) #checked for accuracy
        C=k1-(1+k2)*math.cos(math.radians(i))+k3
        d=(B**2)-(4*A*C) #checked for accuracy
        sqd=math.sqrt(d)
        A2=2*A
        sol1=(-B-sqd)/A2
        sol2=(-B+sqd)/A2
        solved_angle1=math.degrees(math.tanh(sol1))
        #print(solved_angle1)
        solved_angle2=math.degrees(math.tanh(sol2))
        #print(solved_angle2)
        error1=abs(solved_angle1-expected_output_angles[index_counter])
        error2=abs(solved_angle2-expected_output_angles[index_counter])
        if (error1>error2):
            calculated_output_angles[index_counter]=solved_angle2*2
        else:
            calculated_output_angles[index_counter]=solved_angle1*2
        index_counter=index_counter+1
    return calculated_output_angles
def input_angle_generator(Range=5,initial_input_angle=180, final_input_angle=120):
    array_size=abs(int((final_input_angle-initial_input_angle)/Range))
    input_angles=[0.00]*array_size
    index_counter=0
    for i in range(final_input_angle,initial_input_angle, Range):
        input_angles[index_counter]=i
        index_counter=index_counter+1
    return input_angles
        
def error_ploter():
    input_angles=input_angle_generator(Range=5)
    errors_with_3angles=structural_errors()
    errors_with_5angles=structural_errors(numberOfAngles=5)
    plt.plot(input_angles, errors_with_3angles, color='r', label='Error With 3 Angles')
    plt.plot(input_angles, errors_with_5angles, color='g', label='Error With 5 Angles')
    plt.xlabel("Input Angles")
    plt.ylabel("Errors")
    plt.title("Input Angles Against Errors")
    plt.legend()
    plt.show()
    
#input_angles_from_cheb()
#Output_angles_from_cheb()
#calculateKs(5)
#calculate_link_lengths(numberOfAngles=5)
#transmission_angle()
#least_square_method(5)
#calculateKs()
#structural_errors()
error_ploter()


      
        

    
    