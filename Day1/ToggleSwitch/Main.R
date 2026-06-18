### 1. Model Generation

library(epimod)
model.generation(net_fname = "ToggleSwitch.PNPRO")


### 2. Deterministic Analysis
# Initial condition favouring State 1 (A > B)

model.analysis(
  solver_fname = "ToggleSwitch.solver",
  i_time = 0, f_time = 20, s_time = 0.1,
  n_config = 1,
  parameters_fname = "./Input/parameters_state1.csv"
)
display_data(volume = "./")

# Initial condition favouring State 2 (B > A)
model.analysis(
  solver_fname = "ToggleSwitch.solver",
  i_time = 0, f_time = 20, s_time = 0.1,
  n_config = 1,
  parameters_fname = "./Input/parameters_state2.csv"
)
display_data(volume = "./")

### 3. Stochastic Analysis

model.analysis(
  solver_fname  = "ToggleSwitch.solver",
  i_time = 0, f_time = 50, s_time = 0.1,
  n_run  = 50,
  solver_type = "SSA",
  parameters_fname = "./Input/parameters_symmetric.csv"
)
display_data(volume = "./")
stop_display()

