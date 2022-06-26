#!/usr/bin/env clojure

;; A raytracer written in Clojure,
;; based on “Ray Tracing in One Weekend.” at raytracing.github.io

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; SCENE PARAMETERS

(def image-width 400)
(def image-height 225)
(def focal-length 1.0)
(def samples-per-pixel 100)
(def viewport-height 2.0)
(def max-depth 50)

;; UTILITIES

;; mark all args as double typed
(defmacro defn' [name args body]
  (letfn [(f [s] (if (vector? s)
                   (with-meta (mapv f s) (meta s))
                   (with-meta s {:tag 'double})))]
    `(defn ~name ~(f args) ~body)))

(defn' random ^double [a b] (+ a (* (- b a) ^double (rand))))
(defn' clamp ^double [x xmin xmax] (min xmax (max xmin x)))

;; VECTOR MATH

(defn' vec+ [[e0 e1 e2] [v0 v1 v2]] [(+ e0 v0) (+ e1 v1) (+ e2 v2)])
(defn' vec* [[e0 e1 e2] s] [(* e0 s) (* e1 s) (* e2 s)])
(defn' vec*vec [[e0 e1 e2] [v0 v1 v2]] [(* e0 v0) (* e1 v1) (* e2 v2)])

(defn vec-neg [v] (vec* v -1.0))
(defn' vec- [[e0 e1 e2] [v0 v1 v2]] [(- e0 v0) (- e1 v1) (- e2 v2)])

(defn x ^double [v] (v 0)) (defn y ^double [v] (v 1)) (defn z ^double [v] (v 2))

(defn' vec-len-sq ^double [[e0, e1, e2]] (+ (* e0 e0) (* e1 e1) (* e2 e2)))

(defn vec-len ^double [v] (Math/sqrt (vec-len-sq v)))

(defn' vec-lerp [[u0 u1 u2] [v0 v1 v2] t]
  (let [t- (- 1.0 t)]
    [(+ (* t- u0) (* t v0))
     (+ (* t- u1) (* t v1))
     (+ (* t- u2) (* t v2))]))

(defn' cross [[u0 u1 u2] [v0 v1 v2]]
  [(- (* u1 v2) (* u2 v1))
   (- (* u2 v0) (* u0 v2))
   (- (* u0 v1) (* u1 v0))])

;; same as (reduce + (map * u v))
(defn' dot ^double [[u0 u1 u2] [v0 v1 v2]]
  (+ (* u0 v0) (* u1 v1) (* u2 v2)))

(defn vec-unit [v] (vec* v (/ 1.0 (vec-len v)))) ;; normalize to unit length

(defn vec-rand-in-unit-sphere []
  (loop []
    (let [pt [(random -1 +1) (random -1 +1) (random -1 +1)]]
      (if (>= (vec-len-sq pt) 1.0)
        (recur)
        pt))))

(defn vec-rand-in-hemisphere [normal]
  (let [v (vec-rand-in-unit-sphere)]
    (if (pos? (dot v normal))
      v (vec-neg v))))

;; rays are represented by this triple for now
(defn ->ray [origin direction] [:ray origin direction])
(defn ray-origin [[_ origin _]] origin) (defn ray-direction [[_ _ direction]] direction)

(defn ray-at [ray t] (vec+ (ray-origin ray) (vec* (ray-direction ray) t)))

;; MATERIALS

(defn lambertian [albedo]
  (fn [ray pt normal front?]
    (let [scatter-direction (vec+ pt (vec-rand-in-hemisphere normal))
          scattered (->ray pt scatter-direction)]
      [albedo scattered])))

(defn reflect [v n]
  (vec- v (vec* n (* 2.0 (dot v n)))))

(defn metal [albedo]
  (fn [ray pt normal front?]
    (let [reflected (reflect (vec-unit (ray-direction ray)) normal)
          scattered (->ray pt reflected)]
      (when (pos? (dot reflected normal))
        [albedo scattered]))))

(defn refract [uv n etai-over-etat]
  (let [cos-theta  (Math/min (dot (vec-neg uv) n) 1.0)
        r-out-perp (vec* (vec+ uv (vec* n cos-theta)) etai-over-etat)
        r-out-par  (vec* n (- (Math/sqrt (Math/abs (- 1.0 (vec-len-sq r-out-perp))))))]
    (vec+ r-out-perp r-out-par)))

(defn dielectric [refract-idx]
  (fn [ray pt normal front?]
    (let [refraction-ratio (if-not front? refract-idx (/ 1.0 refract-idx))
          unit-dir         (vec-unit (ray-direction ray))
          refracted        (refract unit-dir normal refraction-ratio)
          scattered        (->ray pt refracted)]
      [[1.0 1.0 1.0] scattered])))

;; SCENE OBJECTS

(defmulti hit (fn [m ray t-min t-max] (:type m)))

(defmethod hit :sphere [{:keys [center ^double radius material]} ray t-min t-max]
  (let [oc     (vec- (ray-origin ray) center)
        a      (vec-len-sq (ray-direction ray))
        half-b (dot oc (ray-direction ray))
        c      (- (vec-len-sq oc) (* radius radius))
        discriminant (- (* half-b half-b) (* a c))]
    (when (>= discriminant 0)
      (let [sqrtd (Math/sqrt discriminant)
            root1 (/ (- (- half-b) sqrtd) a)
            root2 (/ (+ (- half-b) sqrtd) a)]
        (when-let [root (cond (<= t-min root1 t-max) root1
                              (<= t-min root2 t-max) root2)]
          (let [p              (ray-at ray root)
                outward-normal (vec* (vec- p center) (/ 1.0 radius))
                front-face     (neg? (dot (ray-direction ray) outward-normal))
                normal         (if front-face outward-normal (vec-neg outward-normal))]
            {:t        root
             :p        p
             :front    front-face
             :normal   normal
             :material material}))))))

(defmethod hit :list [{:keys [elems]} ray t-min t-max]
  (loop [index (dec (count elems))
         found nil
         t-max t-max]
    (if (neg? index)
      found
      (if-let [found (hit (elems index) ray t-min t-max)]
        (recur (dec index) found (:t found))
        (recur (dec index) found t-max)))))

(def origin [0.0 0.0 0.0])
(def aspect-ratio (/ 16.0 9.0))
(def viewport-width (* (double aspect-ratio) (long viewport-height)))
(def horizontal [viewport-width 0.0 0.0])
(def vertical [0.0 viewport-height 0.0])

(defn background [ray]
  (let [unit-direction (vec-unit (ray-direction ray))
        t (* 0.5 (+ (y unit-direction) 1.0))]
    (vec-lerp [1.0 1.0 1.0] [0.5 0.7 1.0] t)))

(defn ray-color [ray hittable ^long depth]
  (if (neg? depth)
    origin
    (if-let [{:keys [^double p normal material front]} (hit hittable ray 0.001 ##Inf)]
      (if-let [[attenuation scattered] (material ray p normal front)]
        (vec*vec attenuation (ray-color scattered hittable (dec depth)))
        origin)
      (background ray))))

(def lower-left-corner
  (-> origin
      (vec- (vec* horizontal 0.5))
      (vec- (vec* vertical 0.5))
      (vec- [0 0 focal-length])))

(defn render-pixel [hittable ^long i ^long j]
  (loop [sample (int samples-per-pixel)
         r 0.0 g 0.0 b 0.0]
    (if (zero? sample)
      [r g b]
      (let [u   (/ (+ (double (rand)) i) (- (long image-width) 1))
            v   (/ (+ (double (rand)) j) (- (long image-height) 1))
            ray (->ray origin (-> lower-left-corner
                                  (vec+ (vec* horizontal u))
                                  (vec+ (vec* vertical v))))
            c   (ray-color ray hittable max-depth)]
        (recur (dec sample) (+ r (x c)) (+ g (y c)) (+ b (z c)))))))

(defn print-ppm [^long width ^long height pixelcallback]
  (println "P3")
  (println width height)
  (println 255)
  (doseq [j (range height)]
    (binding [*out* *err*]
      (println "Rendering row" j))
    (doseq [px (pmap #(pixelcallback % (- height (int j))) (range width))]
      (let [[^double r ^double g ^double b] px
            scale (/ 1.0 (double samples-per-pixel))]
        (println (int (* 256 (clamp (Math/sqrt (* r scale)) 0.0 0.999)))
                 (int (* 256 (clamp (Math/sqrt (* g scale)) 0.0 0.999)))
                 (int (* 256 (clamp (Math/sqrt (* b scale)) 0.0 0.999))))))))

(let [sphere1 {:type :sphere :center [0.0 0.0 -1.0] :radius 0.5 :material (lambertian [0.7 0.3 0.3])}
      sphere2 {:type :sphere :center [0.0 -100.5 -1.0] :radius 100.0 :material (lambertian [0.8 0.7 0.0])}
      sphere3 {:type :sphere :center [-1.0 0.0 -1.0] :radius 0.5 :material (dielectric 1.5) #_(lambertian [0.3 0.7 0.3])}
      sphere4 {:type :sphere :center [+1.0 0.0 -1.0] :radius 0.5 :material (metal [0.8 0.8 0.8])}
      ]
  (def world {:type :list :elems [sphere1 sphere2 sphere3 sphere4]}))

(print-ppm image-width image-height (partial render-pixel world))
(shutdown-agents)

; EOF