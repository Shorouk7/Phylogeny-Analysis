#%%
!pip install tajimas-d

# %%
from tajimas_d import tajimas_d, pi_estimator, watterson_estimator

sequences = ["AAAA", "AAAT", "AAGT", "AAGT"]
theta_tajima = tajimas_d(sequences)
theta_pi = pi_estimator(sequences)
theta_w = watterson_estimator(sequences)

print("Tajima's D:", theta_tajima)
print("Pi estimator:", theta_pi)
print("Watterson estimator:", theta_w)