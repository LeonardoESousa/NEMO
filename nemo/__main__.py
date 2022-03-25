#!/usr/bin/env python3
import sys
import nemo.tools

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
        freqlog = nemo.tools.fetch_file("frequency",['.out', '.log'])
        with open(freqlog, 'r') as f:
            for line in f:
                if 'Entering Gaussian System' in line:
                    gauss = True
                else:
                    gauss = False
                break
        if gauss:             
            print('You are using a Gaussian log file.')
            template = nemo.tools.fetch_file("QChem template",['.out', '.in'])
            rem, cm, spec = nemo.tools.busca_input(template)
        else:    
            rem, cm, spec = nemo.tools.busca_input(freqlog)        
        print('\nThe suggested configurations for you are:\n')
        print(rem)
        change = input('Are you satisfied with these parameters? y or n?\n')
        if change.lower() == 'n':     
            rem2 = ''
            for elem in rem.split('\n'):
                if len(elem.split()) > 1:
                    if '$' not in elem:
                        base = nemo.tools.default(elem, '{} is set to: {}. If ok, Enter. Otherwise, type the correct value. Type del to delete line.\n'.format(elem.split()[0], elem.split()[-1]))
                        if base.lower() == 'del':
                            base = ''
                    else:    
                        base = elem
                    rem2 += base+'\n'
        num_ex = input("How many excited states?\n")
        try:
            num_ex = int(num_ex)
        except:
            nemo.tools.fatal_error("This must be a number! Goodbye!")
        abs_only = input("Prepare input for absorption or fluorescence spectrum only? (y or n)\n")
        if abs_only.lower() == 'y':
            print('Ok, calculations will only be suitable for absorption or fluorescence spectrum simulations!\n')
            header = "$comment\n{}\n$end\n\n$rem\ncis_n_roots             {}\ncis_singlets            true\ncis_triplets            true\ncalc_soc                false\nSTS_MOM                 true".format(spec,num_ex)
        else:
            print('Ok, calculations will be suitable for all spectra and ISC rate estimates!\n')
            header = "$comment\n{}\n$end\n\n$rem\ncis_n_roots             {}\ncis_singlets            true\ncis_triplets            true\ncalc_soc                true\nSTS_MOM                 true".format(spec,num_ex)
        header =  rem.replace('$rem',header)
        header += '$molecule\n{}\n'.format(cm)
        num_geoms = int(input("How many geometries to be sampled?\n"))
        T = float(input("Temperature in Kelvin?\n"))
        if T <= 0:
            nemo.tools.fatal_error("Have you heard about absolute zero? Goodbye!")
        if gauss:
            import lx.tools
            lx.tools.make_ensemble(freqlog, num_geoms, T, header,'$end\n')
        else:    
            nemo.tools.make_ensemble(freqlog, num_geoms, T, header,'$end\n')    
    elif op == '2':
        nemo.tools.batch() 
    elif op == '3':
        nemo.tools.andamento()
    elif op == '4':
        nemo.tools.abort_batch()        
    elif op == '5':
        opc = nemo.tools.detect_sigma()
        tipo = nemo.tools.get_spec()
        nr = nemo.tools.get_nr() 
        print('The spectrum will be run with the following parameters:\n')
        print('Spectrum type: {}'.format(tipo.title()))
        print('Standard deviation of: {:.3f} eV'.format(opc))
        print('Refractive index: {:.3f}\n'.format(nr))
        change = input('Are you satisfied with these parameters? y or n?\n')
        if change.lower() == 'n':
            tipo = input("What kind of spectrum? Type abs (absorption) or emi (emission)\n")
            if tipo != 'abs' and tipo != 'emi':
                nemo.tools.fatal_error('It must be either one. Goodbye!')
            opc = input("What is the standard deviation of the gaussians?\n")
            try:
                opc = float(opc)
            except: 
                nemo.tools.fatal_error("It must be a number. Goodbye!")  
        tipo = tipo[:3]
        if tipo == 'abs':
            estados = nemo.tools.ask_states("Absorption from which state (S0, S1, T1 ..)\n")
            nemo.tools.spectra('abs', estados, nr, opc)
        elif tipo == 'emi':
            estados = nemo.tools.ask_states("Emission from which state (S1, T1, etc)?\n")
            nemo.tools.spectra('emi', estados, nr, opc)    
    elif op == '6':
        state = input('What is the initial state (S1, T1, S2 ...)? Accepts comma separated values Ex: T1,T2\n')
        from nemo.analysis import isc
        states = state.split(',')
        for state in states:
            isc(state)
    elif op == '7':
        from lx.tools import ld
        ld()
    elif op == '8':
        freqlog = nemo.tools.fetch_file("log",['.log','.out'])
        rem, cm, spec = nemo.tools.busca_input(freqlog)
        G, atomos = nemo.tools.pega_geom(freqlog)
        nemo.tools.write_input(atomos,G,'{}\n$molecule\n{}\n'.format(rem,cm),'$end','geom.lx')
        print('Geometry saved in the geom.lx file.')    
    else:
        nemo.tools.fatal_error("It must be one of the options... Goodbye!")


    
if __name__ == "__main__":
    sys.exit(main())        

