# Analysis of the Heli-GPS Data

## Introduction and Data

After lodgepole pines (Pinus contorta) have been attacked by Mountain
Pine Beetle (MPB), they show different stages. The needles of attacked
trees wil turn from green over yellow to red within the first 12 months
(red-attack). It then takes in between 1-3 years until all needles have
fallen to the ground, when the trees remain standing but without needles
(grey-attack) ([Wulder et al.,
2005](#ref-wulderRemoteSensingSurvey2005)). Depending on the type of
stand, according to
([**mitchellFallRateLodgepole199?**](#ref-mitchellFallRateLodgepole199)),
dead trees begin to fall after 3 to 5 years and 80 - 90 % of trees had
fallen after 11 years. This leads to different stages of beetle attack
that will be considered in this study. All green-attack trees detected
in 2021 nad 2020 are now with high probability in the red-attack stage.
Some of the trees that have been detected as red-attack in 2021, will
still classify as red-attack and others as grey-attack. A clear
distinction based on time of detection is not possible without
additional data, since it can take several years before all needles have
fallen of the tree, the GPS locations of MPB-attacked trees provided by
the government of Alberta ([Wulder et al.,
2005](#ref-wulderRemoteSensingSurvey2005)). Collecting the LiDAR data at
the end of May, all green-attack trees that have been attacked in 2019
and before, as well as all red-attack trees detected in 2020 and before,
will now classify as grey attack. Additionally, it can be expected that
most of the trees that have been detected as red-attack in 2011 and
before, now have fallen to the ground ([Schoennagel et al.,
2012](#ref-schoennagelEffectsMountainPine2012)).

Three broad areas of interest have been selected based on availability
of observed MPB attack in the years 2011, 2019, 2020, and 2021, while
also having a dense network of paved roads and gravel roads. As already
shown by various authors before, lodgepole pine (Pinus contorta) is the
most commonly attacked by MPB, thus, areas where chosen, where lodgepole
pine is one of the most dominant species. This information was based on
a visual analysis of raser data containing the most dominant species in
30 m by 30 m grid cells. The government of Alberta conducts annual
Heli-GPS surveys, usually between August 15 and September 15 by flying
helicopters at low hight and mapping the GPS locations of detected
disturbances in therms of red- and green-attack through observers. These
locations have a positional accuracy of =/- 30 m. In most cases only
patches consisting of at least 3 trees are recorded. For the years under
investigation, the minimum number of trees per patch were \_\_\_. For
each of these locations the corresponding cell from the raster data
containing the tree species is determined and with a box filter of 5 by
5 the most common tree species for each point can be determined. This
analysis is performed for each of the three areas of interest and the
results are shown in fig !{fig:freq_species_11_19_20_21}[Frequency of
tree species within the surrounding area of each MPB survey point for
the years 2011, 2019, 2020, and
2021](%22./graphics/freq_species_11_19_20_21.png%22) ![Frequency how
often each tree species was the most dominant within one cell for each
area of interest and
cumulative.](%22./graphics/freq_species_19_20_21_cum.png%22)

Information on accessibility is based on the road network available
through the government of Alberta ([Government of Alberta,
2021a](#ref-GravelRoad2021), [2021b](#ref-PavedRoad2021)), since all
flight locations with the drone need to be accessible from the road for
efficient data collection. As the drone has to stay in line of sight, it
is not possible to fly far away from the road.

# References

<div id="refs" class="references csl-bib-body hanging-indent"
line-spacing="2">

<div id="ref-GravelRoad2021" class="csl-entry">

Government of Alberta. (2021a). *Layer: Gravel Road (20K) (ID: 23)*.
https://geospatial.alberta.ca/titan/rest/services/utility/access/MapServer/23.

</div>

<div id="ref-PavedRoad2021" class="csl-entry">

Government of Alberta. (2021b). *Layer: Paved Road (20K) (ID: 14)*.
https://geospatial.alberta.ca/titan/rest/services/utility/access/MapServer/14.

</div>

<div id="ref-schoennagelEffectsMountainPine2012" class="csl-entry">

Schoennagel, T., Veblen, T. T., Negron, J. F., & Smith, J. M. (2012).
Effects of Mountain Pine Beetle on Fuels and Expected Fire Behavior in
Lodgepole Pine Forests, Colorado, USA. *PLOS ONE*, *7*(1), e30002.
<https://doi.org/10.1371/journal.pone.0030002>

</div>

<div id="ref-wulderRemoteSensingSurvey2005" class="csl-entry">

Wulder, M. A., White, J. C., & Dymond, C. C. (2005). *Remote sensing in
the survey of mountain pine beetle impacts: Review and recommendations*.
Pacific Forestry Centre.

</div>

</div>
