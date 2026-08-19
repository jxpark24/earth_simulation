FROM nginx:1.27-alpine

RUN apk add --no-cache curl

# Vendor three.js and the Earth textures into the image so the site has no
# runtime dependency on any CDN.
WORKDIR /usr/share/nginx/html
RUN mkdir -p vendor assets && \
    curl -sSfL -o vendor/three.module.js \
        https://unpkg.com/three@0.160.0/build/three.module.js && \
    curl -sSfL -o vendor/OrbitControls.js \
        https://unpkg.com/three@0.160.0/examples/jsm/controls/OrbitControls.js && \
    for t in earth_atmos_2048.jpg earth_specular_2048.jpg earth_normal_2048.jpg earth_clouds_1024.png; do \
        curl -sSfL -o "assets/$t" \
          "https://raw.githubusercontent.com/mrdoob/three.js/r160/examples/textures/planets/$t"; \
    done && \
    sed -i "s#from 'three'#from '/vendor/three.module.js'#g" vendor/OrbitControls.js && \
    ls -la vendor assets

COPY nginx.conf /etc/nginx/conf.d/default.conf
COPY index.html .
COPY js/ js/

EXPOSE 80
