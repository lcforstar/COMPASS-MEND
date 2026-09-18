# Predicting Plant Water from Space: Vegetation Optical Depth as a Prognostic Ecosystem Variable

Jizixing Zhu^1^, Xiaodong Gao^1,\*^, Tianjie Zhao^2^, Xining Zhao^1^

^1^ State Key Laboratory of Soil and Water Conservation and Desertification Control, College of Soil and Water Conservation Science and Engineering, Northwest A&F University, Yangling, Shaanxi Province, China

^2^ State Key Laboratory of Remote Sensing Science, Aerospace Information Research Institute, Chinese Academy of Sciences, Beijing, China

^\*^ Corresponding author. E-mail: Xiaodong Gao

**Plant water content governs growth, drought survival, and land–atmosphere carbon and water exchange, and its daily swing carries the fingerprint of stress. Microwave satellites sense it through vegetation optical depth (VOD), but twice-daily passes miss this diurnal signal. We argue the fix is to let process-based models predict VOD as a continuous output, not merely retrieve it from satellites.**

## The observational gap

Vegetation water content (VWC) governs transpiration, drought survival, and wildfire risk, and its diurnal swing is a fingerprint of stress: hydraulic failure is now recognised as a leading cause of tree mortality, and post-drought recovery depends on whether plants can rehydrate (Konings et al., 2019). The need for VWC is long established — but measuring it is the bottleneck. Direct field measurement is destructive, sparse, and point-scale, and plant traits vary too much from place to place to sample stem by stem; crucially, no field campaign delivers a *continuous* record at the scales ecology needs.

This bottleneck is why microwave VOD became the observable of choice. Unlike SIF (photosynthetic activity) or thermal infrared (canopy temperature), VOD responds directly to whole-plant water: it rises with total canopy water — biomass times water content (Wigneron et al., 2007) — covers the globe daily, and at L-band reaches past leaves into branches and stems, where thirst is registered (Momen et al., 2017). VOD thus eases the spatial side of the measurement problem. But it does not solve the temporal side: a satellite passes a given spot only about twice a day, so a VOD retrieval captures a tree in one or two fixed states, not the continuous response — the pre-dawn-to-midday drawdown, the lag behind soil drying — through which stress actually reveals itself (Holtzman et al., 2023). VOD has therefore solved part of the VWC problem and left the most diagnostic part, its time evolution, open.

The way forward, we argue, is to also generate VOD from a process-based model with explicit plant hydraulics. Such a model serves two ends at once: its simulated VWC can be checked against the (sparse) field measurements, while its forward-modelled VOD is continuous in time and can be compared with — or assimilated together with — the satellite signal. Model and observation then reinforce each other, and VOD becomes a continuous ecological observable rather than a twice-daily snapshot (Figure 1).

## Why retrieved VOD is not enough

The way the field uses VOD today is to invert satellite measurements into a retrieved product. As retrieval methods have improved, each advance has exposed a deeper limit. Five stand out.

1. **The VOD–water link is empirical.** Reading VOD as plant water needs a curve that differs by plant type, so the same VOD means different things in different places (Konings et al., 2019). The link is fitted, not derived.

2. **Products disagree on the value.** Sensors differ in frequency, angle, orbit, and quality, so VOD values are mutually biased and must be rescaled before combining, which can create artificial spatial jumps (Moesinger et al., 2024). Two products can give different VOD for the same place, with no physical reference to settle which is right.

3. **A retrieval is a snapshot, not a record.** Modern retrievals reach 9 km (Zhao et al., 2021) and gap-filling yields seamless daily maps (Hu et al., 2025), yet the value is fixed at the overpass moment: it says nothing about how canopy water changes in between, or why.

4. **The orbit limits how often we see it.** Microwave satellites revisit a spot about twice a day, twelve hours apart, while plant water swings between a pre-dawn high and an early-afternoon low about eight hours apart (Holtzman et al., 2023). Two snapshots miss most of that swing — and the swing is exactly where the drought-stress signal lives.

5. **Inversion from VOD to traits is ambiguous.** VOD constrains some plant parameters well and others poorly, so a good fit to one number does not mean the system is identified (Holtzman et al., 2023).

## Why models have not closed the gap either

Process models are the natural place to supply what retrieval cannot. But most do not really compute plant water: they estimate it from leaf area or biomass, or read it from a fixed table, collapsing the flow from soil to stem to leaf into a single scaling constant. Where VOD is generated at all, it is treated as a *diagnostic* — a by-product read off biomass — rather than a *prognostic* state carrying the model's own water physics. Such a model cannot reproduce the delay between soil and canopy drying, or the daily cycle from a well-watered dawn canopy to a thirsty midday one — the very behaviour that distinguishes a stressed plant from a healthy one (Kennedy et al., 2019; L. Li et al., 2021). The shortfall is conceptual: water is present as a factor but never placed at the centre that the satellite signal is checked against.

## The components already exist

The encouraging point is that none of the parts are missing. Models with plant hydraulics already track water through the soil–plant–atmosphere system, linking stomatal conductance, transpiration, xylem vulnerability, and tissue capacitance into a chain that determines the water held in leaves, stems, and roots (Kennedy et al., 2019; L. Li et al., 2021). The output of that chain — total canopy water as biomass times VWC, with diurnal stem refilling and depletion shaping VOD hourly — is exactly what a forward VOD operator needs (Holtzman et al., 2023).

The optical side proves the route works: process models already turn a simulated plant state into the SIF a satellite would see, from SCOPE (van der Tol et al., 2009) and SCOPE 2.0 (Yang et al., 2021) to land-model implementations at site and global scale (R. Li et al., 2022; Miyauchi et al., 2025). SIF and VOD work the same way — both are satellite signals computed forward from a plant state — so what works for SIF can work for VOD. On the microwave side this is no longer only an idea either: recent systems compute VOD from a modelled state to constrain plant hydraulics at a single site (Holtzman et al., 2023), and assimilate VOD globally through an empirical observation operator driven by biomass (Knorr et al., 2025).

A second class of microwave forward operator, detailed discrete-scattering models, takes a different route. These describe the canopy branch by branch and are accurate in principle, but the mechanism is placed in canopy electromagnetic geometry: they demand structural detail — size, angle, and arrangement of scatterers — rarely known at the relevant scales (Ferrazzoli & Guerriero, 1996), and the empirical fits used as substitutes are most fragile precisely under the extreme conditions we most need to characterise. A plant-hydraulics approach places the mechanism instead in the water, energy, and carbon balances the model already solves, so VOD follows from the simulated water through a comparatively simple relation, without the demanding structural inputs.

What remains is to unite these strands — global coverage, prognostic plant water, and a physical forward operator — with water, not biomass, as the driver. This middle ground, between Holtzman's prognostic-but-single-site demonstration and Knorr's global-but-empirical operator, is where the opportunity lies.

## Outlook

Treating VOD as something a model predicts — driven by its own plant water rather than a biomass proxy — opens a clear ecological path. Matching modelled plant water against the satellite signal sharpens what is otherwise unobservable: root water uptake, post-drought recovery, and proximity to hydraulic failure (Knorr et al., 2025). Because the same water must simultaneously satisfy the model's water, energy, and carbon balances, physics replaces the fitted VOD–water link, the model fills the gaps between passes by how water actually moves, and inversion ambiguity shrinks. More frequent passes from paired or geostationary orbits would capture the daily cycle that two-a-day orbits miss, and a model turns those extra views into ecological inference rather than repeated snapshots (Holtzman et al., 2023). Anchoring still needs field measurements of plant water and water potential — in particular paired predawn and midday observations that resolve the diurnal swing the model must reproduce. Done well, prognostic VOD could provide a globally consistent benchmark for ecosystem drought resilience and the stability of the terrestrial carbon sink. Plant water is a state ecology needs to watch; prognostic VOD is how we can watch it from space.

## Figure 1: From a retrieved VOD to a prognostic plant water variable

![Figure 1](figures/vod_paradigm_figure1.png)

## Code and data availability

No new code or data were generated for this perspective. All datasets and models discussed are described in the cited publications and available from their respective public archives.

## Author contributions

JZ conceived the perspective and wrote the manuscript. TZ contributed expertise on microwave VOD retrieval and revised the manuscript. XG supervised the work and revised the manuscript. All authors approved the final version.

## Competing interests

The authors declare no competing interests.

## Acknowledgements

[To be completed: funding sources and grant numbers.]

---

## References

1. Holtzman, N., Wang, Y., Wood, J. D., Frankenberg, C. & Konings, A. G. *Water Resour. Res.* **59**, e2023WR035481 (2023).
2. Konings, A. G., Rao, K. & Steele-Dunne, S. C. *New Phytol.* **223**, 1166–1172 (2019).
3. Kennedy, D. et al. *J. Adv. Model. Earth Syst.* **11**, 485–513 (2019).
4. Li, L. et al. *J. Adv. Model. Earth Syst.* **13**, e2020MS002214 (2021).
5. Wigneron, J.-P. et al. *Remote Sens. Environ.* **107**, 639–655 (2007).
6. Momen, M. et al. *J. Geophys. Res. Biogeosci.* **122**, 3031–3046 (2017).
7. Zhao, T. et al. *Remote Sens. Environ.* **257**, 112321 (2021).
8. Moesinger, L. et al. *Earth Syst. Sci. Data* **16**, 4573–4617 (2024).
9. Hu, D. et al. *Earth Syst. Sci. Data* **17**, 2849–2872 (2025).
10. Ferrazzoli, P. & Guerriero, L. *IEEE Trans. Geosci. Remote Sens.* **34**, 433–443 (1996).
11. van der Tol, C., Verhoef, W., Timmermans, J., Verhoef, A. & Su, Z. *Biogeosciences* **6**, 3109–3129 (2009).
12. Yang, P., Prikaziuk, E., Verhoef, W. & van der Tol, C. *Geosci. Model Dev.* **14**, 4697–4712 (2021).
13. Li, R. et al. *J. Adv. Model. Earth Syst.* **14**, e2021MS002747 (2022).
14. Knorr, W. et al. *Geosci. Model Dev.* **18**, 2137–2159 (2025).
15. Miyauchi, T. et al. *Geosci. Model Dev.* **18**, 2329–2347 (2025).
