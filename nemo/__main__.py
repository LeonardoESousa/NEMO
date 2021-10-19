#!/usr/bin/env python3
import sys
from nemo.tools import *


def main():
    print("#     # ####### #     # #######")
    print("##    # #       ##   ## #     #")
    print("# #   # #       # # # # #     #")
    print("#  #  # #####   #  #  # #     #")
    print("#   # # #       #     # #     #")
    print("#    ## #       #     # #     #")
    print("#     # ####### #     # #######")
    print("----------Photophysics---------\n")
    print("Choose your option:\n")
    print("ENSEMBLE SETUP:")
    print("\t1 - Generate the inputs for the nuclear ensemble calculation")
    print("\t2 - Run the ensemble calculations")
    print("\t3 - Check the progress of the calculations")
    print("\t4 - Abort my calculations")
    print('SPECTRUM SIMULATIONS:')
    print("\t5 - Generate the spectrum")
    print("INTERSYSTEM CROSSING (ISC):")
    print("\t6 - Estimate ISC rates")
    print('EXCITON ANALYSIS:')
    print("\t7 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths")
    print('OTHER FEATURES:')
    print("\t8 - Retrieve last geometry from log file") 
    op = input()
    if op == '1':
        freqlog = fetch_file("frequency",['.out'])
        rem, cm, spec = busca_input(freqlog)        
        print('\nThe suggested configurations for you are:\n')
        print(rem)
        change = input('Are you satisfied with these parameters? y or n?\n')
        if change.lower() == 'n':     
            rem2 = ''
            for elem in rem.split('\n'):
                if len(elem) > 0:
                    if '$' not in elem:
                        base = default(elem, '{} is set to: {}. If ok, Enter. Otherwise, type the correct value\n'.format(elem.split()[0], elem.split()[-1]))
                    else:
                        base = elem
                    rem2 += base+'\n'
        num_ex = input("How many excited states?\n")
        try:
            num_ex = int(num_ex)
        except:
            fatal_error("This must be a number! Goodbye!")
        header = "$comment\n{}\n$end\n\n$rem\ncis_n_roots             {}\ncis_singlets            true\ncis_triplets            true\nRPA                     false\ncalc_soc                true\nSTS_MOM                 true".format(spec,num_ex)
        header =  rem.replace('$rem',header)
        header += '$molecule\n{}\n'.format(cm)
        num_geoms = int(input("How many geometries to be sampled?\n"))
        T = float(input("Temperature in Kelvin?\n"))
        if T <= 0:
            fatal_error("Have you heard about absolute zero? Goodbye!")
        sample_geom(freqlog, num_geoms, T, header)    
    elif op == '2':
        batch() 
    elif op == '3':
        andamento()
    elif op == '4':
        abort_batch()        
    elif op == '5':
        opc = detect_sigma()
        tipo = get_spec()
        nr = get_nr() 
        print('The spectrum will be run with the following parameters:\n')
        print('Spectrum type: {}'.format(tipo.title()))
        print('Standard deviation of: {:.3f} eV'.format(opc))
        print('Refractive index: {:.3f}\n'.format(nr))
        change = input('Are you satisfied with these parameters? y or n?\n')
        if change.lower() == 'n':
            opc = input("What is the standard deviation of the gaussians?\n")
            try:
                opc = float(opc)
            except: 
                fatal_error("It must be a number. Goodbye!")  
            tipo = input("What kind of spectrum? Type abs (absorption) or flu (fluorescence) or pho (phosphorescence)\n")
            if tipo != 'abs' and tipo != 'flu' and tipo != 'pho':
                fatal_error('It must be either one. Goodbye!')
        tipo = tipo[:3]
        gather_data(opc)
        if tipo == 'abs':
            estados = ask_states("How many excited states in the absorption spectrum?\n")
            spectra('abs', estados, nr)
        elif tipo == 'flu':
            estados = ask_states("Fluorescence from which excited state (1,2, etc)?\n")
            spectra('fluor', estados, nr)
        elif tipo == 'pho':    
            estados = ask_states("\nPhosphorescence from which excited state (1,2, etc)?\n")
            spectra('phosph', estados, nr)    
    elif op == '6':
        state = input('What is the initial state? (S1, T1, S2 ...)\n')
        from nemo.analysis import isc
        isc(state)
    elif op == '7':
        from lx.tools import ld
        ld()
    elif op == '8':
        freqlog = fetch_file("log",['.log','.out'])
        rem, cm, spec = busca_input(freqlog)
        G, atomos = pega_geom(freqlog)
        write_input(atomos,G,'{}\n$molecule\n{}\n'.format(rem,cm),'$end','geom.lx')
        print('Geometry saved in the geom.lx file.')    
    else:
        fatal_error("It must be one of the options... Goodbye!")


    
if __name__ == "__main__":
    sys.exit(main())        

