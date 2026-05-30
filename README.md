# Pitman–Yor Process Clustering for Amplification Timing

This project estimates copy-number amplification timing from mutation VAFs and clusters timing events using a truncated Pitman–Yor process mixture model.

---

## Segment Timing

For segment $s$, let

$$C_s = n_{A,s} + n_{B,s}$$

be the total copy number, and let $\rho$ be tumor purity.

The expected VAF for a mutation present on one copy is

$$p_{1s} = \frac{\rho}{C_s \rho + 2(1-\rho)}.$$

The expected VAF for a mutation duplicated by the amplification is

$$p_{2s} = \frac{2\rho}{C_s \rho + 2(1-\rho)}.$$

Given the mean observed VAF

$$\bar{v}_s = \frac{1}{M_s} \sum_{m=1}^{M_s} \mathrm{VAF}_{sm},$$

we estimate

$$\bar{v}_s = \theta_{1s} p_{1s} + \theta_{2s} p_{2s},$$

with

$$\theta_{1s} + \theta_{2s} = 1.$$

Therefore,

$$\theta_{2s} = \frac{\bar{v}_s - p_{1s}}{p_{2s} - p_{1s}}, \qquad \theta_{1s} = 1 - \theta_{2s}.$$

The timing estimate is $\tau_s \in [0,1]$.

For CNLOH and tetrasomy,

$$\tau_s = \frac{2\theta_{2s}}{2\theta_{2s} + \theta_{1s}}.$$

For trisomy,

$$\tau_s = \frac{3\theta_{2s}}{\theta_{1s} + 2\theta_{2s}}.$$

For general amplifications,

$$\tau_s = \frac{C_s \theta_{2s}}{\theta_{1s} + C_s \theta_{2s}}.$$

---

## Pitman–Yor Mixture Model

The segment timings are modeled as

$$\tau_s \mid z_s = k,\, \mu_k,\, \sigma^2 \sim \mathcal{N}(\mu_k, \sigma^2).$$

The marginal density is

$$p(\tau_s) = \sum_{k=1}^{K} w_k \, \mathcal{N}(\tau_s \mid \mu_k, \sigma^2).$$

The random mixing measure is

$$G = \sum_{k=1}^{K} w_k \delta_{\mu_k}.$$

The stick-breaking weights are

$$w_k = v_k \prod_{\ell < k}(1 - v_\ell),$$

with

$$v_k \sim \mathrm{Beta}(1 - d,\, \alpha).$$

Here $\alpha > 0$ is the concentration parameter and $0 \le d < 1$ is the Pitman–Yor discount. When $d = 0$ the model behaves like a Dirichlet-process mixture; when $d > 0$ it supports heavier-tailed cluster-size distributions.

---

## Allocation Probability

For each segment,

$$\Pr(z_s = k \mid {-}) \propto w_k \, \mathcal{N}(\tau_s \mid \mu_k, \sigma^2).$$

After normalization,

$$\Pr(z_s = k \mid {-}) = \frac{w_k \, \mathcal{N}(\tau_s \mid \mu_k, \sigma^2)}{\displaystyle\sum_{h=1}^{K} w_h \, \mathcal{N}(\tau_s \mid \mu_h, \sigma^2)}.$$

---

## Gibbs Updates

Let $n_k = |\{s : z_s = k\}|$ and $n_{>k} = \sum_{h=k+1}^{K} n_h$.

The stick-breaking update is

$$v_k \mid {-} \sim \mathrm{Beta}\left(1 - d + n_k,\; \alpha + dk + n_{>k}\right).$$

For occupied clusters,

$$\mu_k \mid {-} \sim \mathcal{N}\left(\bar{\tau}_k,\; \frac{\sigma^2}{n_k}\right),$$

where

$$\bar{\tau}_k = \frac{1}{n_k} \sum_{s:\, z_s=k} \tau_s.$$

---

## Retrospective Likelihood

Given timing $\tau_s$, recover $\theta_{2s}$.

For CNLOH and tetrasomy,

$$\theta_{2s} = \frac{\tau_s}{2 - \tau_s}.$$

For trisomy,

$$\theta_{2s} = \frac{\tau_s}{3 - 2\tau_s}.$$

For general amplifications,

$$\theta_{2s} = \frac{\tau_s}{C_s - \tau_s(C_s - 2)}.$$

Then $\theta_{1s} = 1 - \theta_{2s}$.

With depth $D = 100$ and alternate read count $a_{sm} = \mathrm{round}(D \cdot \mathrm{VAF}_{sm})$, the mutation likelihood is

$$L_{sm} = \theta_{1s} \, \mathrm{Binomial}(a_{sm} \mid D, p_{1s}) + \theta_{2s} \, \mathrm{Binomial}(a_{sm} \mid D, p_{2s}).$$

The total log-likelihood is

$$\mathcal{L} = \sum_s \sum_m \log(L_{sm} + \varepsilon).$$

---

## Posterior Predictive Check

For each mutation,

$$c_{sm} \sim \mathrm{Categorical}(\theta_{1s}, \theta_{2s}).$$

Then

$$p_{sm} = \begin{cases} p_{1s}, & c_{sm} = 1, \\ p_{2s}, & c_{sm} = 2. \end{cases}$$

Replicated read counts are sampled as

$$a_{sm}^{\mathrm{rep}} \sim \mathrm{Binomial}(D, p_{sm}),$$

with replicated VAF

$$\mathrm{VAF}_{sm}^{\mathrm{rep}} = \frac{a_{sm}^{\mathrm{rep}}}{D}.$$

Observed and replicated VAFs are compared using the Kolmogorov–Smirnov statistic per cluster:

$$D_k = \sup_x \left| F_k^{\mathrm{obs}}(x) - F_k^{\mathrm{rep}}(x) \right|.$$

---

## Model Comparison

The PYP model is compared against TickTack using

$$\Delta\mathcal{L} = \mathcal{L}_{\mathrm{PYP}} - \mathcal{L}_{\mathrm{TickTack}}.$$

The within-cluster variance is

$$\mathrm{WCV} = \frac{\displaystyle\sum_{k=1}^{K} (n_k - 1)\,\widehat{\mathrm{Var}}\{\tau_s : z_s = k\}}{S - K}.$$

The relative improvement is

$$100 \cdot \frac{\mathrm{WCV}_{\mathrm{TickTack}} - \mathrm{WCV}_{\mathrm{PYP}}}{\mathrm{WCV}_{\mathrm{TickTack}}}.$$

