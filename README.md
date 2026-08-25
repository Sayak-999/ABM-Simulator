# Agent-Based Simulator for the COVID-19 Pandemic

An interactive **agent-based model (ABM)** that simulates how an infection like COVID-19 spreads through a small, synthetic community. Instead of tracking population-level averages (as classical SIR/SEIR compartmental models do), this simulator tracks **every individual ("agent")** day by day — each with its own age, sex, occupation, economic status, mobility, immunity, and comorbidity — and lets you watch an outbreak unfolds under different interventions.

The model is delivered as an [R Shiny](https://shiny.posit.co/) web app (`app.R`). You can run it in your browser without writing any code.

The methodology and results are described in the paper **["An AI-Enabled Agent-Based Simulation Platform for Studying COVID-19 Pandemic"](https://arxiv.org/abs/2106.11070)** (arXiv:2106.11070).

> **▶️ Live app: [https://anikburman.shinyapps.io/ABM_simulator/](https://anikburman.shinyapps.io/ABM_simulator/)**
>
> **📄 Paper: [arXiv:2106.11070](https://arxiv.org/abs/2106.11070)**

### 📺 Tutorial video

▶️ **[Watch the tutorial video](app_tutorial.mp4)** — a short walkthrough of how to set up and run a simulation. Click the link to open the video in GitHub's built-in player (with audio).

> The same video is also available inside the app's **Overview** tab.

---

## Table of contents

- [What the simulator does](#what-the-simulator-does)
- [How the model works](#how-the-model-works)
- [Quick start (use the live app)](#quick-start-use-the-live-app)
- [Step-by-step usage guide](#step-by-step-usage-guide)
- [Parameter reference](#parameter-reference)
- [Reading the output plots](#reading-the-output-plots)
- [Running the app locally](#running-the-app-locally)
- [Repository structure](#repository-structure)
- [Known issues](#known-issues)
- [Citation](#citation)
- [Authors](#authors)
- [License](#license)

---

## What the ABM-simulator does

Understanding the dynamics of an outbreak is essential for designing effective control measures. This tool lets you experiment with those dynamics directly. You set up a virtual community of several thousand people, drop one or more initial infections into it, optionally impose lockdowns, and then run the clock forward to see how many people become infected, recover, or die over time.

Because every agent behaves according to its own attributes, the model naturally captures realistic, non-linear effects that are hard to express with equations — for example, the fact that an outbreak starting inside a dense slum cluster behaves very differently from one starting in a sparse area, or that lower testing sensitivity quietly drives up the death toll.

You can use it to explore questions such as:

- How much does a lockdown flatten the curve, and what happens when it is lifted?
- Does it matter *where* the first cases appear — inside a crowded cluster or out in the open?
- How do more (or fewer) asymptomatic cases change the total number of infections and deaths?
- What is the cost, in lives, of poorer diagnostic sensitivity?

## How the model works

Each agent is assigned demographic attributes from **2011 India Census** distributions (age, sex, working status, economic status), plus an immunity probability that depends on age and sex. Agents move around a square map each simulated day; movement patterns depend on age, gender, work status, economic status, and whether the agent is abiding by lockdown.

Infection spreads through **proximity**: every infected agent carries a circular "neighbourhood" (the *radius of infection*). When a susceptible agent enters that neighbourhood, they may become infected with a probability that depends on the infector's disease stage and the susceptible person's circumstances.

An infected agent then progresses through a series of disease states over roughly a four-week course:

`Susceptible → Pre-symptomatic → Early-symptomatic → Late-symptomatic → Recovered / Deceased`

or, for a configurable fraction of cases, an **Asymptomatic** path that still transmits the disease. Whether a late-symptomatic agent survives depends on a mortality risk that rises with age and comorbidity, and on whether they were correctly identified by testing (governed by the *sensitivity* parameter). Recovered agents are assumed immune for the rest of the run, and deceased agents are removed from the population.

## Quick start (use the live app)

1. Open the **[live app](https://anikburman.shinyapps.io/ABM_simulator/)**.
2. Go to the **Run Simulation** tab.
3. Set **"Let's get started!!"** to **Yes**.
4. On the map grid, click to place at least one **Infection Cluster** (where the outbreak begins). This is required — without it, nothing happens.
5. Press **START** and watch the progress bar. When it finishes, the result plots animate to completion.

That's the 30-second version. The full walkthrough is below.

## Step-by-step usage guide

**1. Read the Overview tab.**
The **Overview** tab explains the model and every parameter, and includes a short **tutorial video**. Start here if it's your first time.

**2. Switch to the Run Simulation tab and enable the app.**
Set the **"Let's get started!! (To start click 'Yes')"** radio button to **Yes**. The interactive map grid and the **Population Cluster** / **Infection Cluster** buttons will appear.

**3. Choose your parameters in the left sidebar.**
Set population size, simulation length, proportion of asymptomatics, transmission probability, radius of infection, area spread, and sensitivity. See the [parameter reference](#parameter-reference) for what each one means and sensible ranges.

**4. (Optional) Configure lockdowns.**
Set **"Do You Want Lockdown Options?"** to **Yes**. Enter one or more lockdown **start days** as a comma-separated list (e.g. `20,60`) and a **common lockdown length** in days. The app warns you in red if a lockdown would start or end after the simulation ends.

**5. Place clusters on the map.**
   - Click a point on the grid, then click **Infection Cluster** to mark where the first cases enter. **You need at least one** (two or more gives more meaningful results); you can place up to five.
   - Optionally, click a point and then **Population Cluster** to create a dense pocket (e.g. a slum) that holds ~20% of the population. You can place up to five. Population outside clusters is spread evenly across the map.

**6. Run it.**
Make sure the start toggle is still on **Yes**, then press **START**. Watch the progress bar at the bottom right — a full run takes a little while because every agent is tracked on every day.

**7. Explore the results.**
When the run completes you'll see a confirmation message and four plots: an animated map of agents colored by disease state, plus three time-series charts (overall progression, symptomatic stages, and symptomatic vs. asymptomatic). See [Reading the output plots](#reading-the-output-plots).

**8. Start over.**
Press **RESET** to clear everything and configure a new scenario.

## Parameter reference

| Parameter | What it controls | Suggested range |
|---|---|---|
| **Population Size** | Number of agents in the community. | 10,000–20,000 for realistic results (max 20,000). |
| **Number of Simulations** | Length of the run, in **days**. | up to 200 (default 150). |
| **Proportion of Asymptomatic People** | Fraction of infections that follow the asymptomatic (but still infectious) path. Held constant for the whole run. | 0.01–0.99 (e.g. 0.10, 0.25, 0.50). |
| **Asymptomatic Transmission Probability** | Chance an asymptomatic agent infects a close contact. | 0.01–0.99. |
| **Radius of Infection** | Size of the circular neighbourhood around an infected agent within which transmission can occur. | 0.005–0.25 units. |
| **Area Spread** | Side length of the square region the population lives in (controls density). | 5–10 (a 10×10 area with 10,000 people is the reference scenario). |
| **Sensitivity** | Probability that an infected person is correctly identified by testing. Lower values → more undetected cases and more deaths. | 0.05–1.0 (e.g. 0.75, 0.85, 0.95). |
| **Lockdown Initiation Time(s)** | Comma-separated start days for one or more lockdowns (e.g. `20,60`). | Any day(s) before the run ends. |
| **Common Length of Lockdowns** | Duration in days, applied equally to every lockdown phase. | 5–25 days. |
| **Population Cluster(s)** | Dense pockets (slums) placed by clicking the map. Together they hold ~20% of the population. | 0–5 clusters. |
| **Infection Cluster(s)** | Epicentre(s) where the first cases appear, placed by clicking the map. **At least one required.** | 1–5 (2+ recommended). |

## Reading the output plots

After a run you get four linked visualizations:

- **Agent map** — every agent plotted at its location and colored by disease state (Susceptible, Pre-symptomatic, Early-symptomatic, Late-symptomatic, Asymptomatic, Recovered, Deceased). This is the spatial view of the outbreak.
- **Overall Infection Progression** — daily **Active cases**, cumulative **Recovered**, and cumulative **Deceased** over time. Lockdown windows are shaded and labelled.
- **Infection stages for Symptomatic patients** — how many agents are in the Pre-, Early-, and Late-symptomatic stages on each day.
- **Symptomatic vs. Asymptomatic** — the daily split of active infections between the two paths.

Typical things you'll observe (and which the accompanying paper documents): without lockdown, active cases in a 10,000-person community peak around week 5–6 and fade by about week 15; outbreaks seeded inside dense clusters peak earlier and linger longer; lockdowns flatten and delay the peak but cases can rebound when restrictions lift; lower sensitivity raises deaths; and a higher asymptomatic proportion lowers total cases and deaths.

## Running the app locally

You need **R** (and optionally **RStudio**).

```r
# 1. Install the required packages (once)
install.packages(c(
  "shiny", "shinyjs", "shinyWidgets",
  "ggplot2", "mvtnorm", "shinythemes"
))

# 2. From the repo root, launch the app
shiny::runApp("app.R")
```

The app reads two data files that **must stay in the repo root** next to `app.R`:

- `Urban_age_dist.txt` — urban age distribution (West Bengal, Census 2011)
- `age__sex_im0.txt` — age/sex immunity inputs

The tutorial video lives in the **`www/`** folder and is served automatically by Shiny.

> This project includes a [`packrat/`](https://rstudio.github.io/packrat/) folder capturing the original package versions, so the environment can be restored reproducibly if you prefer.

## Repository structure

```
.
├── app.R                  # The Shiny app: model functions + UI + server
├── Urban_age_dist.txt     # Age distribution input (Census 2011)
├── age__sex_im0.txt       # Age/sex immunity input
├── www/
│   └── app_tutorial.mp4   # In-app tutorial video  (see Known issues)
├── ABM_Covid19_tex/       # Manuscript (LaTeX source, figures, bibliography)
├── packrat/               # Captured R package environment
└── README.md
```

## Known issues

- **Tutorial video filename mismatch.** `app.R` loads the tutorial video as `www/simvid_updated.mp4` (see the `output$tb` block), but the file currently in `www/` is named `app_tutorial.mp4`. Until these match, the video won't display in the Overview tab. Fix by either renaming the file to `simvid_updated.mp4` **or** updating the `src` in `app.R` to `app_tutorial.mp4`.

## Citation

If you use this simulator or model in your work, please cite the arXiv preprint:

> Burman, A., Chatterjee, S., Ghosh, P., & Mukhopadhyay, I. (2021). *A Flexible Agent-Based Model to Study COVID-19 Outbreak — A Generic Approach.* arXiv:2106.11070 [q-bio.PE]. https://doi.org/10.48550/arXiv.2106.11070

A BibTeX entry:

```bibtex
@misc{burman2021flexibleabm,
  title         = {An AI-Enabled Agent-Based Simulation Platform for Studying COVID-19 Pandemic},
  author        = {Burman, Anik and Chatterjee, Sayak and Ghosh, Pramit and Mukhopadhyay, Indranil},
  year          = {2021},
  eprint        = {2106.11070},
  archivePrefix = {arXiv},
  primaryClass  = {q-bio.PE},
  doi           = {10.48550/arXiv.2106.11070},
  url           = {https://arxiv.org/abs/2106.11070}
}
```

📄 **Read the paper:** [arXiv:2106.11070](https://arxiv.org/abs/2106.11070) · [PDF](https://arxiv.org/pdf/2106.11070)

## Authors

- **Anik Burman** — Johns Hopkins University, Department of Biostatistics
- **Sayak Chatterjee** — University of Pennsylvania, Department of Statistics and Data Science
- **Pramit Ghosh** — ICMR-National Institute for Research in Bacterial Infections, Department of Epidemiology
- **Indranil Mukhopadhyay** — University of Nebraska-Lincoln, Department of Statistics

_Anik Burman and Sayak Chatterjee contributed equally to this work._

## License

A license has not been chosen yet. Until a `LICENSE` file is added to this repository, all rights are reserved and reuse is by permission of the authors. _(Consider adding an OSI-approved license such as MIT or GPL-3.0 before publishing.)_

---

> **Setup checklist before pushing:** swap any `<your-username>/<repo-name>` placeholders, add a `LICENSE` file, resolve the tutorial-video filename note above, and (after the first push) drag-drop `app_tutorial.mp4` into the README on github.com and paste the resulting URL over `PASTE_GITHUB_USER_ATTACHMENTS_VIDEO_URL_HERE`.
