import streamlit as st
from astro_math import *
from utils import chunks
from functools import partial
from math import ceil

if not hasattr(st, "solvable_equations"):
    st.solvable_equations = []

if not hasattr(st, "last_solve"):
    st.last_solve = None

if not hasattr(st, "unit_info"):
    st.unit_info = []

if not hasattr(st, "vals_saves"):
    st.vals_saves = {}

st.set_page_config(page_title="AstroSolver")
c1, c2, c3, c4, c5, c6 = st.columns([1 for _ in range(6)])
cols = [c1, c2, c3, c4, c5, c6]
names = [i for i in vals]
names = chunks(names, ceil(len(names)/6))

def get_sol():
    st.solvable_equations = er.solvable(vals)

def key_0(k):
    st.session_state[k] = ""

def all_0():
    for i in vals:
        key_0(i)
    for i in vals:
        set_val(i, None)
    get_sol()

def lamb3(g):
    s = st.session_state[g]
    try:
        s = float(s)
    except:
        pass
    _ = set_val(g, s if s != "" else None) if pre_screen(g, s) or type(s) in [float, int] or s == "" else key_0(g)
    synchronize_from_val(g)
    get_sol()

def lamb4(eq):
    try:
        val, key = eq.solve(vals)
        s = str(val[-1])
        vals[key] = s
        set_val(key, s)
        st.last_solve = (s, measure_of[key])
        set_input_value(*st.last_solve)
        synchronize_from_val(key)
        get_sol()
    except IndexError:
        st.write("No Possible Solutions")

def set_input_value(text, measure):
    st.session_state.input_value = text
    try:
        st.session_state.measure = measure.name
    except Exception as err:
        print(type(err), err)
        st.session_state.measure = "LENGTH"

def set_unit_info():
    try:
        st.unit_info = as_all_units("", ur.parse_expression(st.session_state.input_value), meas=eval("Measure." + st.session_state.measure))
    except DimensionalityError:
        st.write("Invalid Conversion")

def clear_unit_info():
    st.unit_info = []

def synchronize_from_val(val):
    if vals[val] is not None:
        try:
            s = float(st.session_state[val])
            if units[val] is not None:
                st.session_state[val] = "{:~P}".format((vals[val] * units[val]))
            else:
                st.session_state[val] = str(float(vals[val]))
        except:
            if st.session_state[val] != "":
                s = ""
                v = (vals[val] * units[val]).to(ur.parse_expression(st.session_state[val]).units)
                pv = pint_val(v)
                if int(pv) == pv - 0.9999999999999999:
                    pv = int(pv) + 1
                st.session_state[val] = str(pv) + " {:~P}".format(v.units)
            else:
                s = ""
                v = vals[val] * (units[val] if units[val] is not None else 1)
                if type(v) in [float, int]:
                    st.session_state[val] = str(v)
                elif type(v) in [str]:
                    vals[val] = float(vals[val])
                    st.session_state[val] = str(v)
                else:
                    pv = pint_val(v)
                    if int(pv) == pv - 0.9999999999999999:
                        pv = int(pv) + 1
                    st.session_state[val] = str(pv) + " {:~P}".format(v.units)
    else:
        st.session_state[val] = ""

def synchronize_from_vals():
    for i in vals:
        synchronize_from_val(i)

def save():
    st.vals_saves[st.session_state["save_name_val"]] = vals

def load():
    _ = vals
    try:
        vals = st.vals_saves[st.session_state["save_name_val"]]
        synchronize_from_vals()
    except Exception as err:
        print(type(err), err)

for (n, i) in zip(names, cols):
    with i:
        for g in n:
            st.text_input(g, on_change=lamb3, key=g, args=(g,))

st.button("Clear", on_click=all_0)
st.button("test", on_click=lambda: print(vals))
save_button, save_name, load_button = st.columns([1, 7, 1])
with save_button:
    st.button("Save", on_click=save)
with save_name:
    st.text_input("Save Key", key="save_name_val")
with load_button:
    st.button("Load", on_click=load)

st.markdown("_" * 20)

with st.expander("Converter"):
    c__1, c__2, c__3 = st.columns([3, 2, 1])
    with c__1:
        st.text_input("Input Value", key="input_value")
        if st.last_solve is not None:
            st.button("Set to Last Solve", on_click=set_input_value, args=st.last_solve)
    with c__2:
        st.selectbox("Measure Of", [i.name for i in Measure], key="measure")
    with c__3:
        st.button("Convert", on_click=set_unit_info)
        st.button("Clear", on_click=clear_unit_info, key="Clear2")
    for v in st.unit_info:
        pv = pint_val(v)
        if int(pv) == pv - 0.9999999999999999:
            pv = int(pv) + 1
        if len(str(pv)) > 6:
            pv = float("%.6g" % pv)
        st.markdown("*" + str(pv) + "*" + " **{:P}**".format(v.units))

st.markdown("_" * 20)

_, center, _ = st.columns([1, 1, 1])
with center:
    st.button("Get Solvable Equations", on_click=get_sol)

if len(st.solvable_equations) > 0:
    for i in st.solvable_equations:
        with st.form(key=i.latex):
            c1, c2 = st.columns([5, 1])
            with c1:
                st.write(i.name)
                st.latex(i.latex)
            with c2:
                st.write("")
                st.form_submit_button(label="Solve for", on_click=lamb4, args=(i,))
                st.latex(latex(i.get_unknowns(vals)[0]))