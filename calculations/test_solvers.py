from fenics import *


def print_dict(d):
    for k, v in d.items():
        print(f'{k}: {v}')

prms = dict(parameters)
print_dict(prms)

print('\nkrylov_solver')
prms_krylov = dict(prms['krylov_solver'])
print_dict(prms_krylov)

print('\nlu_solver')
prms_lu = dict(prms['lu_solver'])
print_dict(prms_lu)

print('\nform_compiler')
prms_form = dict(prms['form_compiler'])
print_dict(prms_form)