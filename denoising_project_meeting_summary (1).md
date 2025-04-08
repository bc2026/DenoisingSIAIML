
# 🧠 Meeting Summary: Denoising Project – April 4, 2025

## 🎯 Objective
Discuss progress on particle detection and tracking using machine learning, image augmentation, and improving model performance.

---

## 🔑 Key Takeaways
- **Model Performance**: Achieved 97% accuracy in identifying points of interest using CNNs; still working on generalizing this to multiple data points and tracking across frames.
- **Challenge in Tracking**: High point identification accuracy doesn’t guarantee consistent tracking across frames. Contrast and signal variability make consistent particle tracking challenging.
- **Image Augmentation**: Image augmentation is being explored to generate artificial data, helping improve model robustness and provide more training data.
- **Contrast & Signal Intensity**: Varying signal intensity (above, at, or below threshold) can affect particle detection. It’s crucial to retain and learn from contrast information.
- **Biological Context**: Particles (with DNA) interact with proteins on the cell membrane. These interactions affect movement and contrast, which should be accounted for in the model.
- **Separation of Concerns**: Need to decouple particle detection and particle tracking processes for clarity and effectiveness.

---

## 🧪 Technical Insights
- **Model Improvements**: Added layers and adjusted dropout (keo) to optimize performance.
- **Grayscale Conversion**: All frames are converted to grayscale before feeding into the model.
- **Trajectory Validation**: Current code tracks particles but consistency across frames (trajectory) needs validation.
- **Orientation Matters**: Particle orientation (not just movement) can provide insight into interactions with cellular proteins.
- **Dynamic vs. Static Fields**: Augmentation can help differentiate dynamic signals from static background.

---

## ✅ Action Items

| Who | Task | Status |
|-----|------|--------|
| 🧑‍💻 You | Finalize and present image augmentation results within 2–3 days | ⏳ In progress |
| 🧑‍🏫 Mentor | Review shared drive items and validate consistency of tracking | ✅ Ongoing |
| 🧑‍💻 You | Explore online real-time tracking tools (e.g., TensorFlow-based) | ⏳ In progress |
| Team | Use free GPU resources (e.g., Lightning AI, National Computing Association) for model training | ⏳ Researching |
| Team | Separate particle detection and tracking logic in modeling pipeline | ⏳ Discussing approach |
| 🧑‍💻 You | Merge new point-detection model with existing lab tracking code | ⏳ To do |

---

## 💡 Opportunities
- Leverage national GPU resources or services like Lightning AI for faster training.
- Improve model interpretability by explicitly separating particle detection vs. tracking.
- Use orientation and intensity dynamics as additional features in model training.
