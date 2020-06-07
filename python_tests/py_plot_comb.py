import matplotlib 
import matplotlib.pyplot as plt 
ω = [1.360933971037598,] 
σ = [1.35676494395742e-05,]
N = [0.017026369327038975,]
σ_red = [1.0739898028356683e-05,]

plt.xlabel('ω eV')
plt.ylabel('σ TPA cross')
plt.legend()
plt.plot(ω,σ,'-Dk',label='Full')
plt.plot(ω,σ_red,'-or',label='Reduced')
plt.show()