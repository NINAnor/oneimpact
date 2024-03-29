library(usethis)

#----------
# build package
usethis::use_citation()

#----------
# build pkgdown website

# Run once to configure package to use pkgdown
usethis::use_pkgdown()

# Run to build the website
pkgdown::build_site()

# publish
usethis::use_pkgdown_github_pages()

# setup CI R-CMD-check
usethis::use_github_action_check_standard()

#-----------
# Update articles
pkgdown::build_articles()

