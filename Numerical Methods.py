import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from sympy import symbols, sympify, lambdify, diff, solveset, S
import pandas as pd

def bisection_method(f, a, b, tol=1e-6, max_iter=100):
    data = []
    if f(a) * f(b) >= 0:
        return None, data, "f(a) and f(b) must have opposite signs."
    for i in range(1, max_iter+1):
        c = (a + b) / 2
        data.append([i, a, b, c, f(c)])
        if abs(f(c)) < tol or abs(b - a) < tol:
            return c, data, ""
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return c, data, ""

def newton_raphson_method(f, df, x0, tol=1e-6, max_iter=100):
    data = []
    x = x0
    for i in range(1, max_iter+1):
        fx = f(x)
        dfx = df(x)
        if dfx == 0:
            return None, data, "Zero derivative. No solution found."
        x_new = x - fx / dfx
        data.append([i, x, fx, dfx, x_new])
        if abs(x_new - x) < tol:
            return x_new, data, ""
        x = x_new
    return x, data, ""

def secant_method(f, x0, x1, tol=1e-6, max_iter=100):
    data = []
    for i in range(1, max_iter+1):
        if f(x1) - f(x0) == 0:
            return None, data, "Zero denominator. No solution found."
        x2 = x1 - f(x1)*(x1 - x0)/(f(x1) - f(x0))
        data.append([i, x0, x1, x2, f(x2)])
        if abs(x2 - x1) < tol:
            return x2, data, ""
        x0, x1 = x1, x2
    return x2, data, ""

def false_position_method(f, a, b, tol=1e-6, max_iter=100):
    data = []
    if f(a) * f(b) >= 0:
        return None, data, "f(a) and f(b) must have opposite signs."
    for i in range(1, max_iter+1):
        c = (a*f(b) - b*f(a)) / (f(b) - f(a))
        data.append([i, a, b, c, f(c)])
        if abs(f(c)) < tol or abs(b - a) < tol:
            return c, data, ""
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return c, data, ""

def find_multiple_roots(expr):
    x = symbols('x')
    try:
        roots = list(solveset(expr, x, domain=S.Reals))
        roots = [float(r.evalf()) for r in roots if r.is_real]
        return roots
    except Exception:
        return []

def graphical_method(f, x_start=-10, x_end=10, step=0.01):
    x_vals = np.arange(x_start, x_end, step)
    y_vals = f(x_vals)
    roots = []
    for i in range(1, len(x_vals)):
        if np.sign(y_vals[i-1]) != np.sign(y_vals[i]):
            # Linear interpolation for better root estimate
            x0, x1 = x_vals[i-1], x_vals[i]
            y0, y1 = y_vals[i-1], y_vals[i]
            root = x0 - y0 * (x1 - x0) / (y1 - y0)
            roots.append(root)
    return roots, x_vals, y_vals

def incremental_search_method(f, x_start, x_end, dx):
    x_vals = np.arange(x_start, x_end, dx)
    intervals = []
    for i in range(1, len(x_vals)):
        if f(x_vals[i-1]) * f(x_vals[i]) < 0:
            intervals.append((x_vals[i-1], x_vals[i]))
    return intervals

# Update plot_graph to optionally accept x_vals and y_vals for graphical method
def plot_graph(expr, roots, method_name, canvas_frame, x_vals=None, y_vals=None):
    x = symbols('x')
    f = lambdify(x, expr, 'numpy')
    if x_vals is None or y_vals is None:
        x_vals = np.linspace(-10, 10, 400)
        y_vals = f(x_vals)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x_vals, y_vals, label='f(x)')
    ax.axhline(0, color='black', linewidth=0.5)
    for idx, root in enumerate(roots):
        ax.plot(root, f(root), 'ro')
        ax.annotate(f'Root {idx+1}: {root:.5f}', (root, f(root)), textcoords="offset points", xytext=(0,10), ha='center')
    ax.set_title(f'Graphical Representation ({method_name})')
    ax.set_xlabel('x')
    ax.set_ylabel('f(x)')
    ax.legend(['f(x)', 'Roots'])
    ax.grid(True)
    for widget in canvas_frame.winfo_children():
        widget.destroy()
    canvas = FigureCanvasTkAgg(fig, master=canvas_frame)
    canvas.draw()
    canvas.get_tk_widget().pack()
    plt.close(fig)

class NumericalMethodsApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Numerical Methods Root Finder")
        self.create_widgets()

    def create_widgets(self):
        instructions = (
            "Instructions:\n"
            "1. Enter your equation using 'x' as the variable (e.g., x**3 - 2*x - 5).\n"
            "2. Choose a method and enter the required values.\n"
            "3. Click 'Solve' to see the table and graph.\n"
            "Note: Use '**' for powers, '*' for multiplication."
        )
        tk.Label(self.root, text=instructions, justify="left").pack(pady=5)

        frame = tk.Frame(self.root)
        frame.pack(pady=5)

        tk.Label(frame, text="Equation f(x):").grid(row=0, column=0, sticky="e")
        self.equation_entry = tk.Entry(frame, width=40)
        self.equation_entry.grid(row=0, column=1, columnspan=3, sticky="w")

        tk.Label(frame, text="Method:").grid(row=1, column=0, sticky="e")
        self.method_var = tk.StringVar()
        self.method_combo = ttk.Combobox(
            frame,
            textvariable=self.method_var,
            state="readonly",
            values=[
                "Graphical",
                "Incremental Search",
                "Bisection",
                "False Position",
                "Newton-Raphson",
                "Secant"
            ]
        )
        self.method_combo.grid(row=1, column=1, sticky="w")
        self.method_combo.bind("<<ComboboxSelected>>", self.update_inputs)

        # Description label for method info
        self.method_info = tk.Label(self.root, text="", justify="left", wraplength=600, fg="blue")
        self.method_info.pack(pady=5)

        # Input fields for method parameters
        self.param_labels = []
        self.param_entries = []
        for i in range(2):
            lbl = tk.Label(frame, text="")
            ent = tk.Entry(frame, width=10)
            lbl.grid(row=2+i, column=0, sticky="e")
            ent.grid(row=2+i, column=1, sticky="w")
            self.param_labels.append(lbl)
            self.param_entries.append(ent)

        tk.Button(self.root, text="Solve", command=self.solve).pack(pady=5)

        self.table_frame = tk.Frame(self.root)
        self.table_frame.pack(pady=5)

        self.canvas_frame = tk.Frame(self.root)
        self.canvas_frame.pack(pady=5)

    def update_inputs(self, event=None):
        method = self.method_var.get()
        for lbl, ent in zip(self.param_labels, self.param_entries):
            lbl.config(text="")
            ent.delete(0, tk.END)
            ent.config(state="normal")
        if method in ["Bisection", "False Position"]:
            self.param_labels[0].config(text="Lower bound (a):")
            self.param_labels[1].config(text="Upper bound (b):")
        elif method == "Newton-Raphson":
            self.param_labels[0].config(text="Initial guess (x0):")
            self.param_labels[1].config(text="")
            self.param_entries[1].config(state="disabled")
        elif method == "Secant":
            self.param_labels[0].config(text="First guess (x0):")
            self.param_labels[1].config(text="Second guess (x1):")
        elif method == "Incremental Search":
            self.param_labels[0].config(text="Start x:")
            self.param_labels[1].config(text="End x / dx:")
        elif method == "Graphical":
            self.param_labels[0].config(text="")
            self.param_labels[1].config(text="")
            self.param_entries[0].config(state="disabled")
            self.param_entries[1].config(state="disabled")
        else:
            for ent in self.param_entries:
                ent.config(state="disabled")

        # Show method info
        self.method_info.config(text=self.get_method_info(method))

    def get_method_info(self, method):
        info = {
            "Graphical": (
                "What it is for:\n"
                "  - To visually estimate the root(s) of an equation by plotting the function.\n"
                "How it is used:\n"
                "  - The function is plotted over a range of x-values.\n"
                "  - Where the curve crosses the x-axis (f(x) = 0) is an approximate root.\n"
                "  - Useful for getting a rough idea of root locations."
            ),
            "Incremental Search": (
                "What it is for:\n"
                "  - To find intervals where a root exists by checking for sign changes in f(x) over small increments.\n"
                "How it is used:\n"
                "  - Start at an initial x-value and increment by Δx.\n"
                "  - For each interval [x, x+Δx], check if f(x) and f(x+Δx) have opposite signs.\n"
                "  - If they do, a root exists in that interval."
            ),
            "Bisection": (
                "What it is for:\n"
                "  - To accurately find a root within a given interval where f(x) changes sign.\n"
                "How it is used:\n"
                "  - Start with interval [a, b] where f(a) and f(b) have opposite signs.\n"
                "  - Compute midpoint c = (a+b)/2.\n"
                "  - Check the sign of f(c) and update the interval.\n"
                "  - Repeat until the interval is sufficiently small."
            ),
            "False Position": (
                "What it is for:\n"
                "  - To find a root in an interval [a, b] using a linear interpolation, often converges faster than bisection.\n"
                "How it is used:\n"
                "  - Start with interval [a, b] where f(a) and f(b) have opposite signs.\n"
                "  - Compute c using the formula: c = (a*f(b) - b*f(a)) / (f(b) - f(a)).\n"
                "  - Check the sign of f(c) and update the interval.\n"
                "  - Repeat until convergence."
            ),
            "Newton-Raphson": (
                "What it is for:\n"
                "  - To quickly find a root using tangents, requires the derivative of the function.\n"
                "How it is used:\n"
                "  - Start with an initial guess x₀.\n"
                "  - Compute next approximation: xₙ₊₁ = xₙ - f(xₙ)/f'(xₙ).\n"
                "  - Repeat until the value converges."
            ),
            "Secant": (
                "What it is for:\n"
                "  - To find a root like Newton-Raphson, but does not require the derivative.\n"
                "How it is used:\n"
                "  - Start with two initial guesses x₀ and x₁.\n"
                "  - Compute next approximation: xₙ₊₁ = xₙ - f(xₙ)*(xₙ - xₙ₋₁)/(f(xₙ) - f(xₙ₋₁)).\n"
                "  - Repeat until convergence."
            ),
        }
        return info.get(method, "")

    def solve(self):
        eq_str = self.equation_entry.get()
        method = self.method_var.get()
        x = symbols('x')
        try:
            expr = sympify(eq_str)
            f = lambdify(x, expr, 'numpy')
            df = lambdify(x, diff(expr, x), 'numpy')
        except Exception as e:
            messagebox.showerror("Error", "Invalid equation.")
            return

        # Get parameters
        try:
            if method in ["Bisection", "False Position"]:
                a = float(self.param_entries[0].get())
                b = float(self.param_entries[1].get())
            elif method == "Newton-Raphson":
                x0 = float(self.param_entries[0].get())
            elif method == "Secant":
                x0 = float(self.param_entries[0].get())
                x1 = float(self.param_entries[1].get())
            elif method == "Incremental Search":
                x_start = float(self.param_entries[0].get())
                dx = float(self.param_entries[1].get())
                x_end = x_start + 10*dx  # default range if not specified
        except Exception:
            messagebox.showerror("Error", "Invalid input for method parameters.")
            return

        root = None
        data = []
        error_msg = ""
        columns = []
        roots = []
        x_vals = None
        y_vals = None

        if method == "Bisection":
            root, raw_data, error_msg = bisection_method(f, a, b)
            columns = [
                "Iteration", "x_l", "x_u", "x_r", "f(x_l)", "f(x_u)", "f(x_r)", "f(x_l)*f(x_r)", "Remark"
            ]
            data = []
            for i, x_l, x_u, x_r, f_xr in raw_data:
                f_xl = f(x_l)
                f_xu = f(x_u)
                prod = f_xl * f_xr
                if prod > 0:
                    remark = "Go to next interval"
                elif prod < 0:
                    remark = "Revert to x_l, halve"
                else:
                    remark = "Root found"
                data.append([
                    i, round(x_l, 6), round(x_u, 6), round(x_r, 6),
                    round(f_xl, 6), round(f_xu, 6), round(f_xr, 6),
                    round(prod, 6), remark
                ])
        elif method == "Incremental Search":
            x_start = float(self.param_entries[0].get())
            dx = float(self.param_entries[1].get())
            x_end = x_start + 10*dx
            x_vals_arr = np.arange(x_start, x_end, dx)
            columns = [
                "Iteration", "x", "Δx", "f(x)", "f(x+Δx)", "f(x)*f(x+Δx)", "Remark"
            ]
            data = []
            for i in range(1, len(x_vals_arr)):
                x0 = x_vals_arr[i-1]
                x1 = x_vals_arr[i]
                fx0 = f(x0)
                fx1 = f(x1)
                prod = fx0 * fx1
                if prod < 0:
                    remark = f"Root between {round(x0, 6)} and {round(x1, 6)}"
                else:
                    remark = "No sign change"
                data.append([
                    i, round(x0, 6), round(x1-x0, 6), round(fx0, 6), round(fx1, 6), round(prod, 6), remark
                ])
            # For plotting and root display
            intervals = incremental_search_method(f, x_start, x_end, dx)
            roots = [(interval[0] + interval[1]) / 2 for interval in intervals]
            x_vals = x_vals_arr
            y_vals = f(x_vals_arr)
        elif method == "Newton-Raphson":
            root, raw_data, error_msg = newton_raphson_method(f, df, x0)
            columns = [
                "Iteration", "x_n", "f(x_n)", "f'(x_n)", "x_{n+1}", "Remark"
            ]
            data = []
            for i, x_n, fx_n, dfx_n, x_next in raw_data:
                if abs(x_next - x_n) < 1e-6:
                    remark = "Converged"
                elif i == 1:
                    remark = "Initial step"
                else:
                    remark = ""
                data.append([
                    i, round(x_n, 6), round(fx_n, 6), round(dfx_n, 6), round(x_next, 6), remark
                ])
        elif method == "Secant":
            root, raw_data, error_msg = secant_method(f, x0, x1)
            columns = [
                "Iteration", "x_{n-1}", "x_n", "f(x_{n-1})", "f(x_n)", "x_{n+1}", "Remark"
            ]
            data = []
            for i, x_prev, x_curr, x_next, fx_next in raw_data:
                fx_prev = f(x_prev)
                fx_curr = f(x_curr)
                remark = f"Step {i}"
                data.append([
                    i, round(x_prev, 6), round(x_curr, 6), round(fx_prev, 6), round(fx_curr, 6), round(x_next, 6), remark
                ])
        elif method == "False Position":
            root, raw_data, error_msg = false_position_method(f, a, b)
            columns = [
                "Iteration", "x_l", "x_u", "f(x_l)", "f(x_u)", "x_r", "f(x_r)", "Remark"
            ]
            data = []
            for i, x_l, x_u, x_r, f_xr in raw_data:
                f_xl = f(x_l)
                f_xu = f(x_u)
                if abs(f_xr) < 1e-6:
                    remark = "Converged"
                else:
                    remark = "Go to next interval"
                data.append([
                    i, round(x_l, 6), round(x_u, 6), round(f_xl, 6), round(f_xu, 6), round(x_r, 6), round(f_xr, 6), remark
                ])
        elif method == "Graphical":
            roots, x_vals, y_vals = graphical_method(f)
            columns = []
            data = []

        if error_msg:
            messagebox.showerror("Error", error_msg)
            return

        # Show table
        for widget in self.table_frame.winfo_children():
            widget.destroy()
        if data:
            df = pd.DataFrame(data, columns=columns)
            table = ttk.Treeview(self.table_frame, columns=columns, show='headings', height=min(15, len(data)))
            for col in columns:
                table.heading(col, text=col)
                table.column(col, width=110)
            for row in df.values.tolist():
                table.insert('', 'end', values=row)
            table.pack()
        else:
            tk.Label(self.table_frame, text="No iteration data.").pack()

        # Find and display all real roots (symbolically, if possible)
        if method not in ["Graphical", "Incremental Search"]:
            roots = find_multiple_roots(expr)
            if not roots and root is not None:
                roots = [root]

        # Plot graph
        plot_graph(expr, roots, method, self.canvas_frame, x_vals, y_vals)

        # Show roots in a message box
        if roots:
            msg = "\n".join([f"Root {i+1}: {r:.6f}" for i, r in enumerate(roots)])
            messagebox.showinfo("Roots", f"Roots found:\n{msg}")
        else:
            messagebox.showinfo("Roots", "No roots found.")

if __name__ == "__main__":
    root = tk.Tk()
    app = NumericalMethodsApp(root)
    root.mainloop()