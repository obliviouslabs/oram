import sys
import numpy as np
from scipy.stats import linregress

def analyze_tail_bound(log_file):
    counts = {}
    with open(log_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                size = int(parts[0])
                freq = int(parts[1])
                counts[size] = freq
    
    total_samples = sum(counts.values())
    if total_samples == 0:
        print("No valid data found in the log file.")
        return

    max_size = max(counts.keys())
    
    # Calculate P(X >= k)
    # We want to perform linear regression on the tail, where X >= k
    k_vals = []
    log_probs = []
    
    cumulative_count = total_samples
    for k in range(max_size + 1):
        if k not in counts:
            continue
        prob_exceed = cumulative_count / total_samples
        if prob_exceed > 0 and cumulative_count < total_samples: # Ignore the point where P=1 to fit the tail better, and skip 0 probability
            k_vals.append(k)
            log_probs.append(np.log2(prob_exceed))
        
        cumulative_count -= counts[k]
        if cumulative_count <= 0:
            break

    if len(k_vals) < 2:
        print("Not enough data points to perform linear regression.")
        return

    # Perform linear regression: log2(P(X >= k)) = m * k + c
    slope, intercept, r_value, p_value, std_err = linregress(k_vals, log_probs)
    
    print(f"Linear Regression Results:")
    print(f"Slope (m): {slope:.4f}")
    print(f"Intercept (c): {intercept:.4f}")
    print(f"R-squared: {r_value**2:.4f}")
    print(f"Equation: log2(P(X >= k)) = {slope:.4f} * k + {intercept:.4f}")
    print("-" * 40)
    
    # Calculate required stash size for given failure probabilities
    failure_probs_log2 = [-40, -60, -64, -80]
    print("Required Stash Size (k) to achieve specific failure probabilities:")
    for target_log2 in failure_probs_log2:
        # target_log2 = m * k + c
        # k = (target_log2 - c) / m
        required_k = (target_log2 - intercept) / slope
        print(f"Probability 2^{target_log2}: {required_k:.2f} (ceiling: {int(np.ceil(required_k))})")

    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(8, 6))
    plt.plot(k_vals, log_probs, 'bo', label='Experimental Data', alpha=0.7)
    
    # Plot regression line
    k_array = np.array(k_vals)
    fitted_line = slope * k_array + intercept
    plt.plot(k_array, fitted_line, 'r-', label=f'Regression Line ($log_2(P) = {slope:.2f}k {intercept:+.2f}$)')
    
    plt.xlabel('Stash Size (k)')
    plt.ylabel('$\\log_2(P(X \\geq k))$')
    plt.title('Tail Bound Analysis')
    plt.legend()
    plt.grid(True)
    
    output_pdf = log_file.rsplit('.', 1)[0] + '_tail_bound.pdf' if '.' in log_file else log_file + '_tail_bound.pdf'
    plt.savefig(output_pdf)
    print(f"\nGraph saved to {output_pdf}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python tail_bound_analysis.py <logfile>")
        sys.exit(1)
    
    analyze_tail_bound(sys.argv[1])
