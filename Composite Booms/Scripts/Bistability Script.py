import numpy as np
import math
from composipy import OrthotropicMaterial, LaminateProperty, LaminateStrength
import itertools
import json
import os
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except Exception:
    go = None
    make_subplots = None


#Geometry & Material Properties
ply_thickness = 0.125e-3
h             = 2.25e-2 #Flattened Height
w             = 0.5e-2 #Web Height
R             = 21e-3 #Radius of Curvature of Mandrels
alpha         = (h-w)/(2*np.pi*R) #

#Define
class Materials:
    class UD_CF :
        E1 = 164 #GPa
        E2 = 8.55 #GPa
        v12 = 0.274 
        G12 = 3.51 #GPa
        AW  = 56 #gsm
        
    class PW_CF :
       E1  = 76.5 #GPa
       E2  = 76.5 #GPa
       v12 = 0.031
       G12 = 3.79 #GPa
       AW  = 90 #gsm

#Create Materials
UD_Material = OrthotropicMaterial(Materials.UD_CF.E1, Materials.UD_CF.E2, Materials.UD_CF.v12, Materials.UD_CF.G12, ply_thickness)
PW_Material = OrthotropicMaterial(Materials.PW_CF.E1, Materials.PW_CF.E2, Materials.PW_CF.v12, Materials.PW_CF.G12, ply_thickness)





def BistabilityCalculation (StackingArray, Material):
    # Laminate Properties (use passed StackingArray)
    o_stacking = StackingArray[0]
    i_stacking = StackingArray[1]
    a_stacking = StackingArray[2]

    O_Laminate = LaminateProperty(o_stacking, Material)
    I_Laminate = LaminateProperty(i_stacking, Material)
    A_Laminate = LaminateProperty(a_stacking, Material)


    class Laminate_A_Matrices :
        O = O_Laminate.A
        I = I_Laminate.A
        A = A_Laminate.A
    class Laminate_B_Matrices :
        O = O_Laminate.B
        I = I_Laminate.B
        A = A_Laminate.B
    class Laminate_D_Matrices :
        O = O_Laminate.D
        I = I_Laminate.D
        A = A_Laminate.D

    A_Matrices = Laminate_A_Matrices()
    B_Matrices = Laminate_B_Matrices()
    D_Matrices = Laminate_D_Matrices()

    class Laminate_D_star_Matrices :
        O = D_Matrices.O - np.dot(np.matrix.transpose(B_Matrices.O), np.dot( np.linalg.inv(A_Matrices.O) , B_Matrices.O))
        I = D_Matrices.I - np.dot(np.matrix.transpose(B_Matrices.I), np.dot( np.linalg.inv(A_Matrices.I) , B_Matrices.I))
        A = D_Matrices.A - np.dot(np.matrix.transpose(B_Matrices.A), np.dot( np.linalg.inv(A_Matrices.A) , B_Matrices.A))

    D_star_Matrices = Laminate_D_star_Matrices()

    # Non-dimensionalizing
    class ND_D_Matrices :
        D_ref = D_star_Matrices.O[0,0]  + D_star_Matrices.I[0,0] + D_star_Matrices.A[0,0]
        O = D_star_Matrices.O / D_ref
        I = D_star_Matrices.I / D_ref
        A = D_star_Matrices.A / D_ref

    D_cap = ND_D_Matrices()
    b_cap = w / h

    # Bistability Criterion
    C_cap = (h - w)*(D_cap.I[0,1] - D_cap.O[0,1])/(h*(D_cap.O[0,0] + D_cap.I[0,0]) + w*D_cap.A[0,0])

    ddUcap_ddCcap = D_cap.O[0,0]  + D_cap.I[0,0] + w*D_cap.A[0,0]/h

    ddUcap_ddtheta2 = (2*(D_cap.O[0,1] - D_cap.I[0,1])*(h - w)**2) / (h*(h*(D_cap.O[0,0] + D_cap.I[0,0]) + w*D_cap.A[0,0])**2) * (h*((D_cap.O[0,0] + D_cap.I[0,0])*(D_cap.O[1,1] + D_cap.I[1,1]) - 2*D_cap.I[0,1]*(D_cap.I[2,2] + D_cap.O[2,2] + 0.5*D_cap.I[0,1]) + 2*D_cap.O[0,1]*(D_cap.I[0,1] + D_cap.I[2,2] + D_cap.O[2,2] - 0.5*D_cap.O[0,1])) + w*(D_cap.I[0,1]*(D_cap.I[0,1] - 2*D_cap.A[2,2] - 2*D_cap.O[0,1]) + D_cap.O[0,1]*(D_cap.O[0,1] + 2*D_cap.A[2,2]) + D_cap.A[0,0]*(D_cap.I[1,1] + D_cap.O[1,1])))
    bistable = 1 if (ddUcap_ddtheta2 > 0 and C_cap > 0.1) else 0
    r_coil = R/max(C_cap,0.1)
    r_coil_cm = r_coil * 1e2
    
    # Return a compact result dict (convert numpy scalars to Python floats)
    return {
        'bistable': int(bistable),
        'ddUcap_ddtheta2': float(ddUcap_ddtheta2),
        'C_cap': float(C_cap),
        'r_coil': float(r_coil_cm)
    } 
    
Ply_Angles = [-45, 0, 45, 90]

def run_batch(ply_angles, outer_counts=(2), inner_counts=(3), materials=None, a_stack=None, out_file=None):
    if materials is None:
        materials = {'UD': UD_Material}
    if a_stack is None:
        a_stack = [0]
    results = []
    for mname, mat in materials.items():
        for oc in outer_counts:
            for ic in inner_counts:
                # iterate over ordered combinations (per-ply angles)
                for o_stack in itertools.product(ply_angles, repeat=oc):
                    for i_stack in itertools.product(ply_angles, repeat=ic):
                        stacking_array = [list(o_stack), list(i_stack), list(a_stack)]
                        res = BistabilityCalculation(stacking_array, mat)
                        entry = {
                            'material': mname,
                            'outer_count': oc,
                            'inner_count': ic,
                            'outer_stack': list(o_stack),
                            'inner_stack': list(i_stack),
                            'bistable': res['bistable'],
                            'ddUcap_ddtheta2': res['ddUcap_ddtheta2'],
                            'C_cap': res['C_cap'],
                            'r_coil': res['r_coil']
                        }
                        results.append(entry)

    if out_file:
        # ensure directory exists
        out_dir = os.path.dirname(out_file)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)

        # collect only bistable entries
        bistable_entries = [r for r in results if r.get('bistable') == 1]

        # load existing entries if file exists and is a JSON list
        existing = []
        if os.path.exists(out_file):
            try:
                with open(out_file, 'r') as f:
                    existing = json.load(f)
                if not isinstance(existing, list):
                    existing = []
            except Exception:
                existing = []

        # append new bistable entries, avoid exact-duplicate dicts
        for entry in bistable_entries:
            if entry not in existing:
                existing.append(entry)

        with open(out_file, 'w') as f:
            json.dump(existing, f, indent=2)

    return results


def plot_interactive_3d(results, save_html=None):
    """Create an interactive 3D scatter: outer config (x), inner config (y), C_cap (z).
    Results is a list of dict entries produced by `run_batch` or loaded from JSON.
    """
    if go is None:
        print("Plotly not available. Install with: pip install plotly")
        return None

    # Keep only bistable results
    results = [r for r in results if int(r.get('bistable', 0)) == 1]
    if not results:
        print("No bistable results to plot.")
        return None

    # Represent stacks as short strings for labels
    outer_strs = [str(r['outer_stack']) for r in results]
    inner_strs = [str(r['inner_stack']) for r in results]

    # Create categorical-to-numeric mappings
    unique_outer = sorted(list(dict.fromkeys(outer_strs)))
    unique_inner = sorted(list(dict.fromkeys(inner_strs)))
    outer_map = {s: i for i, s in enumerate(unique_outer)}
    inner_map = {s: i for i, s in enumerate(unique_inner)}

    x = [outer_map[s] for s in outer_strs]
    y = [inner_map[s] for s in inner_strs]
    z = [r.get('r_coil', float('nan')) for r in results]

    color = [r.get('ddUcap_ddtheta2', 0.0) for r in results]
    size = [6 if r.get('bistable', 0) else 4 for r in results]

    hover = [
        f"outer={r['outer_stack']}<br>inner={r['inner_stack']}<br>r_coil={r.get('r_coil'):.6f}<br>ddU={r.get('ddUcap_ddtheta2'):.6f}<br>bistable={r.get('bistable')}"
        for r in results
    ]

    scatter = go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='markers',
        marker=dict(
            size=size,
            color=color,
            colorscale='Sunsetdark',
            colorbar=dict(title='ddUcap_ddtheta2'),
            opacity=0.9
        ),
        hovertext=hover,
        hoverinfo='text'
    )

    fig = go.Figure(data=[scatter])
    fig.update_layout(
        scene=dict(
            xaxis=dict(title='Outer stack', tickmode='array', tickvals=list(outer_map.values()), ticktext=list(outer_map.keys())),
            yaxis=dict(title='Inner stack', tickmode='array', tickvals=list(inner_map.values()), ticktext=list(inner_map.keys())),
            zaxis=dict(title='R_coil')
        ),
        margin=dict(l=0, r=0, b=0, t=30),
        title='Bistability: Outer vs Inner vs C_cap'
    )

    if save_html:
        try:
            fig.write_html(save_html)
            print(f"Interactive plot saved to {save_html}")
        except Exception as e:
            print(f"Failed to save plot HTML: {e}")

    return fig


def plot_2d_curves(results, save_html=None):
    """Create two 2D curves (stacked):
    - Top: C_cap vs outer stack (lines grouped by inner stack)
    - Bottom: ddUcap_ddtheta2 vs outer stack (lines grouped by inner stack)
    """
    if go is None or make_subplots is None:
        print("Plotly not available. Install with: pip install plotly")
        return None

    # Keep only bistable results
    results = [r for r in results if int(r.get('bistable', 0)) == 1]
    if not results:
        print("No bistable results to plot (no bistable entries found).")
        return None

    outer_strs = [str(r['outer_stack']) for r in results]
    inner_strs = [str(r['inner_stack']) for r in results]

    unique_outer = sorted(list(dict.fromkeys(outer_strs)))
    unique_inner = sorted(list(dict.fromkeys(inner_strs)))

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.12,
                        subplot_titles=("C_cap vs Outer stack (grouped by Inner stack)",
                                        "ddUcap_ddtheta2 vs Outer stack (grouped by Inner stack)"))

    for inner in unique_inner:
        y_C = []
        y_ddU = []
        for outer in unique_outer:
            matches = [r for r in results if str(r['outer_stack']) == outer and str(r['inner_stack']) == inner]
            if matches:
                # average if multiple
                vals_C = [m.get('C_cap', float('nan')) for m in matches]
                vals_ddU = [m.get('ddUcap_ddtheta2', float('nan')) for m in matches]
                y_C.append(float(np.nanmean(vals_C)))
                y_ddU.append(float(np.nanmean(vals_ddU)))
            else:
                y_C.append(float('nan'))
                y_ddU.append(float('nan'))

        fig.add_trace(go.Scatter(x=unique_outer, y=y_C, mode='lines+markers', name=f"inner={inner}"), row=1, col=1)
        fig.add_trace(go.Scatter(x=unique_outer, y=y_ddU, mode='lines+markers', name=f"inner={inner}"), row=2, col=1)

    fig.update_xaxes(title_text='Outer stack', row=2, col=1)
    fig.update_yaxes(title_text='C_cap', row=1, col=1)
    fig.update_yaxes(title_text='ddUcap_ddtheta2', row=2, col=1)
    fig.update_layout(height=700, showlegend=True, title_text='Bistability 2D curves')

    if save_html:
        try:
            fig.write_html(save_html)
            print(f"2D curves saved to {save_html}")
        except Exception as e:
            print(f"Failed to save 2D curves HTML: {e}")

    return fig






if __name__ == '__main__':
    # Run the full batch and save to JSON next to this script
    script_dir = os.path.dirname(__file__)
    outpath = os.path.join(script_dir, 'bistability_results.json')
    results = run_batch(Ply_Angles, outer_counts=(2,3), inner_counts=(2,3), out_file=outpath)
    bistable_count = sum(r['bistable'] for r in results)
    total = len(results)
    print(f"Computed {total} cases, Bistable count = {bistable_count}. Results saved to {outpath}")
    # Create an interactive 3D plot (outer stack vs inner stack vs C_cap)
    plot_path = os.path.join(script_dir, 'bistability_plot.html')
    try:
        fig = plot_interactive_3d(results, save_html=plot_path)
    except Exception as e:
        print(f"Plotting failed: {e}")
    # Create 2D curve plots and save HTML
    plot2d_path = os.path.join(script_dir, 'bistability_2d.html')
    try:
        fig2 = plot_2d_curves(results, save_html=plot2d_path)
    except Exception as e:
        print(f"2D plotting failed: {e}")