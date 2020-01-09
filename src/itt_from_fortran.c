#include "ittnotify.h"

__itt_domain* itt_domain_create(char *name) {
  return __itt_domain_create(name); 
}

void itt_frame_begin(__itt_domain *domain, __itt_id *id) {
  __itt_frame_begin_v3(domain, id); 
}

void itt_frame_end(__itt_domain *domain, __itt_id *id) {
  __itt_frame_end_v3(domain,id);
}
