function utility = util(c,l)

load params
utility = ((c - l^omega/omega)^(1-sigm)-1)/(1-sigm);