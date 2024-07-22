# Modern Metagames
# Version 0.01
# Author: Teppichporsche

library(tidyverse)
library(gamlss)
library(broom.mixed)
library(patchwork)
library(stats4)
library(jsonlite)

set.seed(3910)
theme_set(theme_bw())

# Data was generated using the parser with these arguments:
# dotnet MTGOArchetypeParser.App.dll json detect format=Modern startdate=2021-07-01 exclude=League filter=Modern  

# Paper Metagame ----
 
# clean data ----

json <- read_json("mtgo_modern.json")
tib_dat <- tibble(results = json$Data)

data_full <- tib_dat |> hoist(results, 
                Tournament = "Tournament", 
                Wins = "Wins", 
                Losses = "Losses",
                Meta = "Meta",
                url = "AnchorUri",
                Archetype = c("Archetype", "Archetype"))

data_full <- data_full |> rename_with(tolower)

# select only tournaments that were not reported through mtgo

data_paper <- data_full |> 
  filter(!grepl("www.mtgo.com", url))

table(data_paper$meta)
table(data_paper$tournament)

# create number of games 
data_paper <- data_paper |> 
  mutate(wins = as.numeric(wins)) |> 
  mutate(losses = as.numeric(losses)) |> 
  filter(!is.na(wins)) |> 
  filter(!is.na(losses)) |> 
  mutate(games = wins + losses) |> 
  filter(!is.na(games))

# Create archetypes:
arch_paper <- data_paper |> 
  filter(archetype != "Unknown") |> 
  group_by(archetype, meta) |> 
  summarise(allwins = sum(wins, na.rm = TRUE),
            alllosses = sum(losses, na.rm = TRUE),
            games = allwins + alllosses) |> 
  mutate(winp = allwins/games) |> 
  filter(!is.na(winp))

# We have 2505 archetypes, most of which have very little games registered. We
# need some regularization.

# First, we can use a global prior to regularize the data

# log-likelihood function
ll <- function(alpha, beta) {
  x <- arch_paper$allwins
  total <- arch_paper$games
  -sum(VGAM::dbetabinom.ab(x, total, alpha, beta, log = TRUE))
}

# maximum likelihood estimation
m <- mle(ll, start = list(alpha = 1, beta = 10), method = "L-BFGS-B",
         lower = c(0.0001, .1))
ab <- coef(m)

alpha0 <- ab[1]
beta0 <- ab[2]

# We roughly get a alpha = 82, beta = 85 distribution
# The result looks like this

arch_paper |> 
  ggplot() +
  geom_density(aes(winp), linewidth = 1) +
  stat_function(fun = function(x) dbeta(x, alpha0, beta0), color = "red",
                linewidth = 1) +
  xlab("Win % average")

# Clearly, our prior misses some important features of the data, notably in the 
# lower win %. But it might still do a good job regularizing.
# However, it would scale up the lower win % decks, which might not be realistic.


# Individual Priors ----

# We can do better than that by giving each archetype its own prior, which is influenced
# by the number of games played. We assume that good archetypes are played more
# often than very bad archetypes.

# Create MLE estimates for each archetype in each metagame.
# We could fit a more complex function instead of the log-linear effect, but this
# seems to overfit in my testing.

# Metas also can be larger or smaller depending on releases, so we fit seperate 
# slopes for each meta.

fit <- gamlss(cbind(allwins, games - allwins) ~ log(games)*meta,
              data = arch_paper,
              family = BB(mu.link = "identity"))
td <- tidy(fit)
td 
# The number of games played clearly has an influence on win% 

## extract fitted values

mu <- fitted(fit, parameter = "mu")
sigma <- fitted(fit, parameter = "sigma")

arch_paper$mu <- mu
arch_paper$sigma <- sigma

arch_paper <- arch_paper |> 
  mutate(alpha0 = mu / sigma,
         beta0 = (1 - mu) / sigma,
         alpha1 = alpha0 + allwins,
         beta1 = beta0 + games - allwins,
         eb_w = alpha1 / (alpha1 + beta1))

# Plotting Results ----

# This plot shows how the new estimate is a shrunken version of the original
# win% estimates

shrinkage <- ggplot(arch_paper, aes(winp, eb_w, color = games)) +
  geom_point(alpha = 0.3) +
  xlab("Original Win%") +
  ylab("EB Estimate w/ AB term") +
  scale_color_continuous(trans = "log", breaks = 10 ^ (0:4)) +
  xlim(0,1) +
  ylim(0.3,0.7)

shrinkage

ggsave(filename = 'shrinkage.png', plot = shrinkage, width = 7, height = 7, dpi = 300)

# Comparing the number of games and the estimated win% 
rawp <- ggplot(arch_paper, aes(games, winp)) +
  geom_point(alpha = 0.2) +
  xlab("Games") +
  ylab("Win %") +
  scale_x_continuous(trans = "log", breaks = 10 ^ (0:4)) 

hierarchicalp <- ggplot(arch_paper, aes(games, eb_w)) +
  geom_point(alpha = 0.2) +
  xlab("Games") +
  ylab("EB Estimate") +
  scale_x_continuous(trans = "log", breaks = 10 ^ (0:4))

comb <- rawp + hierarchicalp
comb
# We can se the strong amount of regularization for many of these parameters here.
# A lot of decks with low number of games end up at 40% winrate now. The number
# of strong decks is drastically reduced.

ggsave(filename = 'comp_plot.png', plot = comb, width = 7, height = 5, dpi = 300)


# Hall of Fame ----

# What decks very impressive in their metagames? Here, we extract decks which were 
# reliably (using the 99% credible interval) above 50% in their respective metagames

arch_paper <- arch_paper |> 
  mutate(low  = qbeta(.005, alpha1, beta1),
         high = qbeta(.995, alpha1, beta1))

hall_of_fame_MH2 <- arch_paper |> 
  filter(meta!="PostModernHorizons3") |> 
  filter(meta!="PostAssassinsCreed") |> 
  filter(low > .5)
hall_of_fame_MH2


allbest_MH2 <- hall_of_fame_MH2 |> 
  mutate(archetype = reorder(archetype, eb_w)) |> 
  ggplot(aes(eb_w, archetype)) +
  geom_pointrange(aes(xmin = low, xmax = high)) +
  geom_vline(xintercept = 0.5, color = "red", lty = 2) +
  scale_x_continuous(breaks = c(0.5,0.55,0.6,0.65,0.7),limits = c(0.5,0.70)) + 
  xlab("Regularized winrate (w/ 99% credible interval)") +
  ylab("Archetypes") +
  facet_wrap(~meta)

allbest_MH2
ggsave(filename = 'allbest_MH2.png', plot = allbest_MH2, width = 8, height = 15, dpi = 300)

hall_of_fame_MH3 <- arch_paper |> 
  filter(meta=="PostModernHorizons3" | meta=="PostAssassinsCreed") |> 
  filter(low > .5)
hall_of_fame_MH3

allbest_MH3 <- hall_of_fame_MH3 |> 
  mutate(archetype = reorder(archetype, eb_w)) |> 
  ggplot(aes(eb_w, archetype)) +
  geom_pointrange(aes(xmin = low, xmax = high)) +
  geom_vline(xintercept = 0.5, color = "red", lty = 2) +
  scale_x_continuous(breaks = c(0.5,0.55,0.6,0.65,0.7),limits = c(0.5,0.70)) + 
  xlab("Regularized winrate (w/ 99% credible interval)") +
  ylab("Archetypes") +
  facet_wrap(~meta)

allbest_MH3
ggsave(filename = 'allbest_MH3.png', plot = allbest_MH3, width = 6, height = 6, dpi = 300)
